
function [recon,objectiveValues] = mri_jointSparseRecon( data, traj, sMaps, ...
  lambda, varargin )
  % recon = mri_lowRankRecon( data, traj, sMaps, lambda [, 
  %   'alpha', alpha, 'W', W, 'nC', nC ] )
  %
  % Inputs: 
  % data - an NxC array of data values where N is the number of data points and
  %   C is the number of coils
  % traj - a complex vector of N length specifying the k-space trajectory
  % sMaps - sensitivity maps
  % lambda - regularization parameter
  %
  % Optional Inputs:
  % alpha - oversampling factor in Gridding
  % W - window width (in pixels)
  % nC - number of points to use in convolution kernel
  %
  % Written by Nicholas Dwork, Copyright 2020

  p = inputParser;
  p.addParameter( 'alpha', [], @ispositive );
  p.addParameter( 'W', [], @ispositive );
  p.addParameter( 'nC', [], @ispositive );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;

  sSMaps = size( sMaps );
  sImg = sSMaps(1:2);
  if numel( sSMaps ) < 3
    nCoils = 1;
  else
    nCoils = sSMaps(3);
  end
  if size( data, 2 ) ~= nCoils, error( 'Wrong input dimensions' ); end
  nTimes = size( data, 3 );
  nTraj = size( traj, 1 );
  if size( data, 1 ) ~= nTraj, error( 'Wrong input dimensions' ); end

  function out = applyGi( x, type )
    if nargin > 1 && strcmp( type, 'transp' )
      out = iGridT( x, traj, sImg, 'alpha', alpha, 'W', W, 'nC', nC );
    elseif nargin < 2 || strcmp( type, 'notransp' )
      out = iGrid( x, traj, 'alpha', alpha, 'W', W, 'nC', nC );
    end
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

  nData = numel( data );
  function out = applyA( x, op )
    if nargin < 2 || strcmp( op, 'notransp' )
      senseImgs = bsxfun( @times, sMaps, img );
      out = zeros( nTraj, nCoils );
      for coilIndx = 1 : nCoils
        coilImg = senseImgs( :, :, coilIndx );
        out(:,coilIndx) = applyGi( coilImg, 'notransp' );
      end

    elseif nargin > 1 && strcmp( op, 'transp' )
      kData = reshape( x, [ nTraj nCoils ] );
      out = zeros( sImg );
      for coilIndx = 1 : nCoils
        tmp = applyGi( kData(:,coilIndx), 'transp' );
        out = out + sMaps(:,:,coilIndx) .* tmp;
      end

    else
      error( 'Unrecognized op' );
    end

    out = out(:) / sqrt( nData );
  end

  function out = applyAWh( in, op )
    if nargin < 1 || strcmp( op, 'notransp' )
      y = reshape( in, sImg );
      Why = wavTransform( y, 'transp' );
      AWhy = applyA( Why );
      out = AWhy(:);

    elseif strcmp( op, 'transp' )
      y = reshape( in, size( data ) );
      Ahy = applyA( y, 'transp' );
      WAhy = wavTransform( Ahy, 'notransp' );
      out = WAhy(:);

    else
      error( 'Unrecognized op' );
    end
  end

  doCheckAdjoint = false;
doCheckAdjoint = true;
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

    tmp = zeros( [ sImg nTimes ] );
    [checkA,adjErrA] = checkAdjoint( tmp, @applyA );
    if checkA ~= true
      error([ 'applyA failed adjoint test with error ', num2str(adjErrA) ]);
    end
    
    tmp = zeros( [ sImg nTimes ] );
    [checkAWh,adjErrAWh] = checkAdjoint( tmp, @applyAWh );
    if checkAWh ~= true
      error([ 'applyAWh failed adjoint test with error ', num2str(adjErrAWh) ]);
    end
  end

  b = data(:) / sqrt( nData );
  function out = g( in )
    Ain = applyAWh( in );
    out = 0.5 * norm( Ain(:) - b ).^2;
  end

  AWhTb = applyAWh( b, 'transp' );
  function out = gGrad( in )
    AWhin = applyAWh( in );
    AWhTAWhin = applyAWh( AWhin, 'transp' );  
    out = AWhTAWhin - AWhTb;
  end

  function out = proxth( in, t )
    X = reshape( in, [ sImg nTimes ] );
    out = proxL2L1( X, t * lambda );
    out = out( : );
  end

  function out = h( in )
    X = reshape( in, [ prod( sImg ) nTimes ] );
    out = normL2L1( X );
  end

  % Setup inputs to fista_wLS
  innerProd = @(x,y) real( dotP( x, y ) );

  initializeWithGridding = false;
  if initializeWithGridding == true
    weights = makePrecompWeights_2D( traj, sImg, 'alpha', alpha, 'W', W, 'nC', nC );

    y0 = cell( [ 1 1 nTimes ] );
    parfor gTimeIndx = 1 : nTimes
      timeData = data(:,:,gTimeIndx);
      coilRecons = mri_gridRecon( timeData, traj, sImg, 'alpha', alpha, ...
        'W', W, 'nC', nC, 'weights', weights );
      thisRecon = sum( coilRecons .* conj( coilRecons ), 3 );
      %x0(:,:,gTimeIndx) = thisRecon;
      y0{1,1,gTimeIndx} = thisRecon;
    end
    y0 = cell2mat( y0 );
  else
    y0 = zeros( [ sImg nTimes ] );
  end

  % minimize || A Wh y - b ||_2^2 + lambda || y ||_{2,1}
  % Equivalently, minimize g(y) + h(y) where
  % g(y) = || A Wh y - b ||_2^2  and  h(y) = lambda || y ||_{2,1}

  [yStar,objectiveValues] = fista_wLS( y0(:), @g, @gGrad, @proxth, 'N', 30, ...
    't0', 1d-2, 'innerProd', innerProd, 'minStep', 1d-11, 'h', @h, 'verbose', true );

  yStar = reshape( yStar, [ sImg nTimes ] );
  recon = wavTransform( yStar, 'transp' );
end
