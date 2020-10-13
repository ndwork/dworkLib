
function [recon,objectiveValues] = mri_lowRankRecon( data, traj, sMaps, ...
  lambda, varargin )
  % recon = mri_lowRankRecon( data, traj, sMaps, lambda [, 
  %   'alpha', alpha, 'W', W, 'nC', nC ] )
  %
  % Dynamic MRI reconstruction using low rank regularization in the temporal dimension
  % Written as problem (6) in "Accelerated Dynamic MRI Exploiting Sparsity and
  %   Low-Rank Structure: k-t SLR" by Goud et al.
  % Also written according to Zhao, Bo, et al. "Low rank matrix recovery for real-time
  %   cardiac MRI." IEEE, 2010.
  %
  % Inputs: 
  % data - an NxCxT array of data values where N is the number of data points, 
  %   C is the number of coils and T is the number of time points
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

  function out = applyK( x, type )
    if nargin > 1 && strcmp( type, 'transp' )
      kData = reshape( x, [ nTraj nCoils ] );
      out = zeros( sImg );
      for coilIndx = 1 : nCoils
        tmp = applyGi( kData(:,coilIndx), 'transp' );
        out = out + sMaps(:,:,coilIndx) .* tmp;
      end

    elseif nargin < 2 || strcmp( type, 'notransp' )
      img = reshape( x, sImg );
      senseImgs = bsxfun( @times, sMaps, img );
      out = zeros( nTraj, nCoils );
      for coilIndx = 1 : nCoils
        coilImg = senseImgs( :, :, coilIndx );
        out(:,coilIndx) = applyGi( coilImg, 'notransp' );
      end
    end

    out = out(:);
  end

  function out = applyA( x, type )
    if nargin > 1 && strcmp( type, 'transp' )
      kData = reshape( x, [ nTraj nCoils nTimes ] );
      out = zeros( [ sImg nTimes ] );
      for timeIndx = 1 : nTimes
        tmp = applyK( kData(:,:,timeIndx), 'transp' );
        out(:,:,timeIndx) = reshape( tmp, sImg );
      end

    elseif nargin < 2 || strcmp( type, 'notransp' )
      imgs = reshape( x, [ sImg nTimes ] );
      out = zeros( nTraj, nCoils, nTimes );
      for timeIndx = 1 : nTimes
        tmp = applyK( imgs(:,:,timeIndx), 'notransp' );
        out(:,:,timeIndx) = reshape( tmp, [ nTraj nCoils ] );
      end
    end

    out = out(:) / sqrt( nData );
  end

  doCheckAdjoint = false;
  if doCheckAdjoint == true
    tmp = zeros( sImg );
    [checkGi,adjErrGi] = checkAdjoint( tmp, @applyGi );
    if checkGi ~= true
      error([ 'applyGi failed adjoint test with error ', num2str(adjErrGi) ]);
    end

    tmp = zeros( sImg );
    [checkK,adjErrK] = checkAdjoint( tmp, @applyK );
    if checkK ~= true
      error([ 'applyK failed adjoint test with error ', num2str(adjErrK) ]);
    end

    tmp = zeros( [ sImg nTimes ] );
    [checkA,adjErrA] = checkAdjoint( tmp, @applyA );
    if checkA ~= true
      error([ 'applyA failed adjoint test with error ', num2str(adjErrA) ]);
    end
  end

  nData = numel( data );
  b = data / sqrt( nData );
  function out = g( in )
    Ain = applyA( in );
    out = 0.5 * norm( Ain(:) - b(:) ).^2;
  end

  ATb = applyA( b, 'transp' );
  function out = gGrad( in )
    Ain = applyA( in );
    ATAin = applyA( Ain, 'transp' );  
    out = ATAin - ATb;
  end

  function out = proxth( in, t )
    X = reshape( in, [ prod( sImg ) nTimes ] );
    out = proxNucNorm( X, t * lambda );
    out = out( : );
  end

  function out = h( in )
    X = reshape( in, [ prod( sImg ) nTimes ] );
    out = nucNorm( X );
  end

  % Setup inputs to fista_wLS
  innerProd = @(x,y) real( dotP( x, y ) );

  % Intialize x0 with gridding result
  initializeWithGridding = false;
  if initializeWithGridding == true
    weights = makePrecompWeights_2D( traj, sImg, 'alpha', alpha, 'W', W, 'nC', nC );
    x0 = cell( [ 1 1 nTimes ] );
    parfor gTimeIndx = 1 : nTimes
      timeData = data(:,:,gTimeIndx);
      coilRecons = mri_gridRecon( timeData, traj, sImg, 'alpha', alpha, ...
        'W', W, 'nC', nC, 'weights', weights );
      thisRecon = sum( coilRecons .* conj( coilRecons ), 3 );
      x0{1,1,gTimeIndx} = thisRecon;
    end
    x0 = cell2mat( x0 );
  else
    x0 = zeros( [ sImg nTimes ] );
  end
  
  [xStar,objectiveValues] = fista_wLS( x0(:), @g, @gGrad, @proxth, 'N', 50, ...
    't0', 1d-2, 'innerProd', innerProd, 'minStep', 1d-11, 'h', @h, 'verbose', true );

  recon = reshape( xStar, [ sImg nTimes ] );
end
