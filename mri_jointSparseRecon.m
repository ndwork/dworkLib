
function [recon,objValues] = mri_jointSparseRecon( data, traj, sImg, lambda, ...
  varargin )
  % recon = mri_nucNormPLoraksRecon( data, traj, lambda, sImg [, ...
  %   'alpha', alpha, 'W', W, 'nC', nC ] )
  %
  % Performs reconstructions while regularizing the joint sparsity
  % of the images from each coil
  %
  % Inputs: 
  % data - an NxC array of data values where N is the number of data points and
  %   C is the number of coils
  % traj - a complex vector of N length specifying the k-space trajectory
  % lambda - regularization parameter
  % sImg - the size of the image to create
  %
  % Optional Inputs:
  % alpha - oversampling factor in Gridding
  % W - window width (in pixels)
  % nC - number of points to use in convolution kernel
  %
  % Written by Nicholas Dwork, Copyright 2020

  p = inputParser;
  p.addParameter( 'alpha', [], @ispositive );
  p.addParameter( 'debug', false, @islogical );
  p.addParameter( 'W', [], @ispositive );
  p.addParameter( 'nC', [], @ispositive );
  p.addParameter( 'verbose', true, @(x) islogical(x) || isnumeric(x) );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  debug = p.Results.debug;
  W = p.Results.W;
  nC = p.Results.nC;
  verbose = p.Results.verbose;

  sData = size( data );
  nCoils = sData( 2 );
  nTraj = size( traj, 1 );
  if sData(1) ~= nTraj, error( 'Wrong input dimensions' ); end

  function out = applyGi( x, op )
    if nargin > 1 && strcmp( op, 'transp' )
      out = iGridT( x, traj, sImg, 'alpha', alpha, 'W', W, 'nC', nC );
    elseif nargin < 2 || strcmp( op, 'notransp' )
      out = iGrid( x, traj, 'alpha', alpha, 'W', W, 'nC', nC );
    end
  end

  nData = numel( data );
  function out = applyA( x, op )
    if nargin < 2 || strcmp( op, 'notransp' )
      out = cell( 1, nTraj );
      for coilIndx = 1 : nCoils
        out{coilIndx} = applyGi( x(:,:,coilIndx), 'notransp' );
      end
      out = cell2mat( out );
    elseif strcmp( op, 'transp' )
      out = cell(1,1,nCoils);
      for coilIndx = 1 : nCoils
        out{coilIndx} = applyGi( x(:,coilIndx), 'transp' );
      end
      out = cell2mat( out );
    end

    out = out / sqrt( nData );
  end

  wavSplit = zeros(4);  wavSplit(1,1) = 1;
  function out = wavTransform( in, op )
    out = zeros( size( in ) );
    if nargin < 2 || strcmp( op, 'notransp' )
      for coilIndx = 1 : nCoils
        out(:,:,coilIndx) = wtDeaubechies2( in(:,:,coilIndx), wavSplit );
      end
    elseif strcmp( op, 'transp' )
      for coilIndx = 1 : nCoils
        out(:,:,coilIndx) = iwtDeaubechies2( in(:,:,coilIndx), wavSplit );
      end
    end
  end

  function out = applyAW( x, op )
    if nargin < 2 || strcmp( op, 'notransp' )
      x = reshape( x, [ sImg nCoils ] );
      Ax = applyA( x );
      Wx = wavTransform( x );
      out = [ Ax(:); Wx(:); ];
      
    elseif strcmp( op, 'transp' )
      x1 = reshape( x( 1 : nData ), size( data ) );
      ATx = applyA( x1, op );

      x2 = reshape( x( nData + 1 : end ), [ sImg nCoils ] );
      Whx = wavTransform( x2, 'transp' );

      out = ATx(:) + Whx(:);

    else
      error( 'Unrecognized op value' );
    end
  end

  doCheckAdjoint = false;
  if doCheckAdjoint == true
    tmp = zeros( sImg );
    [checkGi,adjErrGi] = checkAdjoint( tmp, @applyGi );
    if checkGi ~= true
      error([ 'applyGi failed adjoint test with error ', num2str(adjErrGi) ]);
    end

    tmp = zeros( [ sImg nCoils ] );
    [checkWav,adjErrWav] = checkAdjoint( tmp, @wavTransform );
    if checkWav ~= true
      error([ 'wavTransform failed adjoint test with error ', num2str(adjErrWav) ]);
    end

    tmp = zeros( [ sImg nCoils ] );
    [checkA,adjErrA] = checkAdjoint( tmp, @applyA );
    if checkA ~= true
      error([ 'applyA failed adjoint test with error ', num2str(adjErrA) ]);
    end

    tmp = zeros( [ sImg nCoils ] );
    [checkAW,adjErrAW] = checkAdjoint( tmp, @applyAW );
    if checkAW ~= true
      error([ 'applyAL failed adjoint test with error ', num2str(adjErrAW) ]);
    end
  end

  b = data(:) / sqrt( nData );
  function out = g( in )
    in1 = in( 1 : nData );
    g1 = 0.5 * norm( in1(:) - b ).^2;

    in2 = in( nData + 1 : end );
    in2 = reshape( in2, [ sImg nCoils ] );
    g2 = lambda * normL2L1( in2 );

disp([ '     Cost parts g1 / g2: ', num2str(g1), ' / ', num2str(g2) ]);

    out = g1 + g2;
  end

  initializeWithGridding = false;
  if initializeWithGridding == true
    weights = makePrecompWeights_2D( traj, sImg, 'alpha', alpha, 'W', W, 'nC', nC );
    x0 = mri_gridRecon( data, traj, sImg, 'alpha', alpha, ...
      'W', W, 'nC', nC, 'weights', weights );
  else
    x0 = zeros( [ sImg nCoils ] );
  end

  function out = proxgConj( v, t )
    v1 = v( 1 : nData );
    out1 = proxConjL2Sq( v1, t, b );

    v2 = v( nData + 1 : end );
    v2 = reshape( v2, [ sImg nCoils ] );
    out2 = proxConjL2L1( v2, lambda );

    out = [ out1(:); out2(:); ];
  end

  f = @(x) 0;
  proxf = @(x,t) x;

  % TODO:  Instead of solving for x, I could solve for y = Wx and then use FISTA
  beta = 1;
  if debug == true
    nIterCP = 10;
  else
    nIterCP = 100;
  end
  if nargout > 1
    [xStar,objValues] = chambollePockWLS( x0(:), proxf, @proxgConj, ...
      'A', @applyAW, 'beta', beta, 'N', nIterCP, 'verbose', verbose, ...
      'doCheckAdjoint', doCheckAdjoint, 'f', f, 'g', @g );
  else
    xStar = chambollePockWLS( x0(:), proxf, @proxgConj, ...
      'A', @applyAW, 'beta', beta, 'N', nIterCP, 'verbose', verbose, ...
      'doCheckAdjoint', doCheckAdjoint );
  end

  recon = reshape( xStar, [ sImg nCoils ] );
end
