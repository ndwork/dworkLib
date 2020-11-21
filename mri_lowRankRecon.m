
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

  nData = numel( data );
  function out = applyA( x, type )
    if nargin < 2 || strcmp( type, 'notransp' )
      %imgs = reshape( x, [ sImg nTimes ] );
      out = cell( 1, 1, nTimes );
      parfor timeIndx = 1 : nTimes
        senseImgs = bsxfun( @times, sMaps, x(:,:,timeIndx) );
        tmp = zeros( nTraj, nCoils );
        for coilIndx = 1 : nCoils
          coilImg = senseImgs( :, :, coilIndx );
          tmp(:,coilIndx) = iGrid( coilImg, traj, 'alpha', alpha, 'W', W, 'nC', nC );
        end
        out{ 1, 1, timeIndx } = tmp;
      end
      out = cell2mat( out );

    elseif nargin > 1 && strcmp( type, 'transp' )
      %kData = reshape( x, [ nTraj nCoils nTimes ] );
      out = cell( 1, 1, nTimes );
      parfor timeIndx = 1 : nTimes
        tmp = zeros([ sImg nCoils ]);
        for coilIndx = 1 : nCoils
          tmp(:,:,coilIndx) = iGridT( x(:,coilIndx,timeIndx), traj, sImg, ...
            'alpha', alpha, 'W', W, 'nC', nC );
        end
        out{ 1, 1, timeIndx } = bsxfun( @times, tmp, sMaps );
      end
      out = cell2mat( out );

    end
    out = out / sqrt( nData );
  end

  doCheckAdjoint = false;
  if doCheckAdjoint == true
    tmp = zeros( [ sImg nTimes ] );
    [checkA,adjErrA] = checkAdjoint( tmp, @applyA );
    if checkA ~= true
      error([ 'applyA failed adjoint test with error ', num2str(adjErrA) ]);
    end
  end

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

  nPixels = prod( sImg );
  function out = proxth( in, t )
    X = reshape( in, [ nPixels nTimes ] );
    out = proxNucNorm( X, t * lambda );
    out = reshape( out, [ sImg nTimes ] );
  end

  function out = h( in )
    X = reshape( in, [ nPixels nTimes ] );
    out = lambda * nucNorm( X );
  end

  innerProd = @(x,y) real( dotP( x, y ) );
  x0 = zeros( [ sImg nTimes ] );

  ATA = @(x) applyA( applyA( x ), 'transp' );
  L = powerIteration( ATA, ones( size( x0 ) ), 'symmetric', true );
  minStep = 0.95 / L;
  t0 = 10 / L;

  nIter = 100;
  [xStar,objectiveValues] = fista_wLS( x0, @g, @gGrad, @proxth, 'N', nIter, ...
    't0', t0, 'innerProd', innerProd, 'minStep', minStep, 'h', @h, 'verbose', true );

  recon = reshape( xStar, [ sImg nTimes ] );
end
