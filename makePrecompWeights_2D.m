
function [weights,nOptIter,flag,res] = makePrecompWeights_2D( kTraj, varargin )
  % [weights,nOptIter,flag,res] = makePrecompWeights_2D( kTraj [, 'sImg', sImg, ...
  %   'alpha', alpha, 'W', W, 'nC', nC, 'alg', alg, 'subAlg', subAlg ] )
  %
  % Determine the density pre-compensation weights to be used in gridding
  %
  % Inputs:
  %   kTraj is a Mx2 element array specifying the k-space trajectory.
  %     The first/second column is the ky/kx location.
  %     The units are normalized to [-0.5,0.5)
  %
  % Optional Inputs:
  %   alpha - the oversampling factor > 1
  %   W - the window width in pixels
  %   nC - the number of points to sample the convolution kernel
  %   alg - a string specifying the algorithm to use
  %     CLS - constrained least squares on trajectory points
  %     FP - specifies fixed point iteration
  %     FP_slow - uses less memory but processes much slower
  %     GP - gradient projection algorithm
  %     GP_sparse - sparse matrix approximation of GP algorithm
  %     JACKSON - one iteration of FP method
  %     LSQR - least squares on trajectory points
  %     SAMSANOV - Constrainted Least Squares on grid points
  %     VORONOI (default) - uses the area of each voronoi cell as the metric of density
  %   gamma - distance weighting parameter for Gradient Projection algorithm
  %     ( default is 0.25 * N )
  %   nIter - number of iterations for iterative algorithms
  %   psfMask - only used by space domain optimizations
  %   sImg is a 2 element array [Ny Nx] representing the number of grid points
  %     A scalar can be supplied and the image is assumed square
  %   verbose - true/false
  %
  % Outputs:
  %   weights - 1D array with density compensation weights
  %
  % Optional Outputs:
  %   nOptIter - number of iterations performed by optimization
  %   flag - flag describing results of optimization (see lsqr
  %     documentation)
  %   res - residual
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp( 'Usage: [weights,nOptIter,flag,res] = makePrecompWeights_2D( kTraj, ...' );
    disp( '  ''sImg'', sImg, ''alpha'', alpha, ''W'', W, ''nC'', nC, ''alg'', alg ] )' );
    disp( '  ''subAlg'', subAlg ' );
    if nargout > 0, weights = []; end
    if nargout > 1, nOptIter = []; end
    if nargout > 2, flag = []; end
    if nargout > 3, res = []; end
    return
  end

  if ~isreal( kTraj ), kTraj = [ real( kTraj(:) ) imag( kTraj(:) ) ]; end

  defaultAlg = 'JACKSON';
  p = inputParser;
  p.addParameter( 'alpha', [], @isnumeric );
  p.addParameter( 'gamma', [], @isnumeric );
  p.addParameter( 'mu', 0, @isnumeric );
  p.addParameter( 'nC', [], @isnumeric );
  p.addParameter( 'nIter', [], @ispositive );
  p.addParameter( 'sImg', [], @ispositive );
  p.addParameter( 'W', [], @ispositive );
  p.addParameter( 'alg', defaultAlg, @(x) true );
  p.addParameter( 'subAlg', [], @(x) true );
  p.addParameter( 'verbose', false, @(x) islogical(x) || isnumeric(x) );
  p.parse( varargin{:} );
  alg = p.Results.alg;
  alpha = p.Results.alpha;
  gamma = p.Results.gamma;
  mu = p.Results.mu;
  nC = p.Results.nC;
  nIter = p.Results.nIter;
  sImg = p.Results.sImg;
  subAlg = p.Results.subAlg;
  verbose = p.Results.verbose;
  W = p.Results.W;
  
  if numel( alg ) == 0, alg = defaultAlg; end
  if numel( sImg ) == 1, sImg = [ sImg sImg ]; end
  if numel( subAlg ) == 0  &&  strcmp( alg, 'GP' ), subAlg = 'proxGrad_wExtrap'; end

  if strcmp( alg, 'JACKSON' ) && numel( sImg ) == 0
    error( 'JACKSON algorithm requires sImg' );
  end

  nOptIter = 1;
  flag = 0;
  res = 0;
  switch alg
    case 'CLS'
      % Constrained least squares density compensation
      [weights,nOptIter,flag,res] = makePrecompWeights_2D_CLS( kTraj, sImg );

    case 'FP_slow'
      % Fixed point iteration
      [weights,flag,res] = makePrecompWeights_2D_FP_slow( kTraj, sImg, 'verbose', verbose, ...
        'alpha', alpha, 'W', W, 'nC', nC );

    case 'FP'
      % Fixed point iteration where matrix is created (may take a lot of memory)
      [weights,flag,res] = makePrecompWeights_2D_FP( kTraj, sImg, 'verbose', verbose, ...
        'alpha', alpha, 'W', W, 'nC', nC );

    case 'GP'
      % Gradient Projection method
      [weights,nOptIter,flag,res] = makePrecompWeights_2D_GP( kTraj, sImg, gamma, mu, nIter, ...
        'alg', subAlg );

    case 'GP_slow'
      % Gradient Projection method
      [weights,nOptIter,flag,res] = makePrecompWeights_2D_GP( kTraj, sImg, gamma, mu, nIter, ...
        'alg', subAlg, 'fast', false );

    case 'GP_sparse'
      % sparse approximation to GP method
      [weights,nOptIter,flag,res] = makePrecompWeights_2D_GP_sparse( kTraj, sImg, gamma, mu, nIter );
      
    case 'JACKSON'
      weights = makePrecompWeights_2D_JACKSON( kTraj, sImg, 'alpha', alpha, 'W', W, 'nC', nC );

    case 'LSQR_slow'
      [weights,nOptIter,flag,res] = makePrecompWeights_2D_LSQR_slow( kTraj, sImg, ...
        'alpha', alpha, 'W', W, 'nC', nC );

    case 'LSQR'
      [weights,nOptIter,flag,res] = makePrecompWeights_2D_LSQR( kTraj, sImg, ...
        'alpha', alpha, 'W', W, 'nC', nC );


    case 'SAMSANOV'
      [weights,nOptIter,flag,res] = makePrecompWeights_2D_SAMSANOV( kTraj, sImg, ...
        'alpha', alpha, 'W', W, 'nC', nC );

    case 'VORONOI'
      weights = makePrecompWeights_2D_VORONOI( kTraj );

    otherwise
      error('makePrecompWeights: Algorithm not recognized');
  end

end


function [weights,nIter,flag,residual] = makePrecompWeights_2D_CLS( traj, N, varargin )
  % Optimization analog of FP with a non-negativity constraint

  p = inputParser;
  p.addParameter( 'alpha', [], @isnumeric );
  p.addParameter( 'W', [], @ispositive );
  p.addParameter( 'nC', [], @ispositive );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;
  
  Ny = 2 * N(1);   [kCy,Cy,~] = makeKbKernel( Ny, Ny, 'alpha', alpha, 'W', W, 'nC', nC );
  Nx = 2 * N(2);   [kCx,Cx,~] = makeKbKernel( Nx, Nx, 'alpha', alpha, 'W', W, 'nC', nC );

  CyCx = Cy(:) * Cx(:)';
  dkCy = zeros( size( kCy ) );  dkCy(1) = kCy(2);  dkCy(2:end) = kCy(2:end) - kCy(1:end-1);
  dkCx = zeros( size( kCx ) );  dkCx(1) = kCx(2);  dkCx(2:end) = kCx(2:end) - kCx(1:end-1);
  dkCykCx = dkCy(:) * dkCx(:)';
  scaling = 1 / sum( CyCx(:) .* dkCykCx(:) );
  Cy = Cy / scaling;
  Cx = Cx / scaling;

  function out = applyA( in, op )
    if nargin < 2 || strcmp( op, 'notransp' )
      out = applyC_2D( in, traj, [ Ny Nx ], kCy, kCx, Cy, Cx );
    elseif strcmp( op, 'transp' )
      out = applyCT_2D( in, traj, [ Ny Nx ], kCy, kCx, Cy, Cx );
    else
      error( 'Unrecognized operator for A' );
    end
  end

  w0 = ( 1 / size(traj,1) ) * ones( size(traj,1), 1 );

  doCheckAdjoint = false;
  if doCheckAdjoint == true
    [checkAdjA,errCheckAdjA] = checkAdjoint( w0, @applyA );
    if checkAdjA == false
      error([ 'Check of A Adjoint failed with error: ', num2str( errCheckAdjA ) ]);
    else
      disp( 'Check of A Adjoint passed' );
    end
  end

  h = @(x) indicatorFunction( x, [ 0, Inf ] );
  proxth = @(x,t) max( x, 0 );

  function out = g( in )
    Ain = applyA( in );
    out = 0.5 * norm( Ain(:) - 1 )^2;
  end

  function out = applyATA( in )
    out = applyA( applyA( in ), 'transp' );
  end

  b = ones( [ Ny Nx ] );
  ATb = applyA( b, 'transp' );
  function out = gGrad( in )
    out = applyA( applyA( in ), 'transp' ) - ATb;
  end

  tol = 1d-6;

  normATA = powerIteration( @applyATA, w0, 'maxIters', 100, ...
    'symmetric', true, 'verbose', true );
  t = 10 / normATA;  % can be too large because line search should take care of it
  minStep = 0.999 / normATA;

verbose = true;
  if verbose == true
    [weights,objValues,relErrs] = fista_wLS( w0, @g, @gGrad, proxth, 'h', h, ...
      'minStep', minStep, 't0', t, 'tol', tol, 'verbose', verbose );   %#ok<ASGLU>
%figure; plotnice( relErrs );
%figure;  semilogynice( relErrs );
  else
    weights = fista_wLS( w0, g, gGrad, proxth );
  end

  if nargout > 1, nIter = numel( relErrs ); end
  flag = 0;
  if nargout > 3, residual = g( weights ); end
end


function [weights,flag,res] = makePrecompWeights_2D_FP_slow( traj, N, varargin )
  % Fixed point iteration defined in "Sampling Density Compensation in MRI:
  % Rationale and an Iterative Numerical Solution" by Pipe and Menon, 1999.

  checknum = @(x) numel(x) == 0 || ( isnumeric(x) && isscalar(x) && (x >= 1) );
  p = inputParser;
  p.addParameter( 'alpha', [], checknum );
  p.addParameter( 'W', [], checknum );
  p.addParameter( 'nC', [], checknum );
  p.addParameter( 'verbose', false, @(x) islogical(x) || isnumeric(x) );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;
  verbose = p.Results.verbose;

  nIter = 8;

  % Make the Kaiser Bessel convolution kernel
  Ny = N(1);  Nx = N(2);
  Gy = Ny;
  [kCy,Cy,~] = makeKbKernel( Gy, Ny, 'alpha', alpha, 'W', W, 'nC', nC );
  Gx = Nx;
  [kCx,Cx,~] = makeKbKernel( Gx, Nx, 'alpha', alpha, 'W', W, 'nC', nC );

  weights = makePrecompWeights_2D_JACKSON( traj, N, 'scaleIt', false );
  for iteration = 2 : nIter
    if verbose ~= 0
      disp(['makePrecompWeights_2D_FP working on iteration ', num2str(iteration) ]);
    end

    oldWeights = weights;
    denom = applyC_2D( oldWeights, traj, traj, kCy, kCx, Cy, Cx );
    weights = oldWeights ./ denom;
  end

  flag = 0;
  res = norm( weights - oldWeights, 2 ) / norm( weights, 2 );

  weights = weights ./ sum( weights(:) );
end


function [weights,flag,res] = makePrecompWeights_2D_FP( traj, N, varargin )
  % Fixed point iteration defined in "Sampling Density Compensation in MRI:
  % Rationale and an Iterative Numerical Solution" by Pipe and Menon, 1999.

  checknum = @(x) numel(x) == 0 || ( isnumeric(x) && isscalar(x) && (x >= 1) );
  p = inputParser;
  p.addParameter( 'alpha', [], checknum );
  p.addParameter( 'W', [], checknum );
  p.addParameter( 'nC', [], checknum );
  p.addParameter( 'verbose', false, @(x) islogical(x) || isnumeric(x) );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;
  verbose = p.Results.verbose;

  nIter = 8;

  % Make the Kaiser Bessel convolution kernel
  Ny = N(1);  Nx = N(2);
  Gy = Ny;
  [kCy,Cy,~] = makeKbKernel( Gy, Ny, 'alpha', alpha, 'W', W, 'nC', nC );
  Gx = Nx;
  [kCx,Cx,~] = makeKbKernel( Gx, Nx, 'alpha', alpha, 'W', W, 'nC', nC );

  C = makeC_2D( traj, traj, kCy, kCx, Cy, Cx );

  nTraj = size( traj, 1  );
  weights = ones( nTraj, 1 );
  for iteration = 1 : nIter
    if iteration == nIter
      oldWeights = weights;
    end
    weights = weights ./ ( C * weights );
  end

  flag = 0;
  res = norm( weights - oldWeights, 2 ) / norm( weights, 2 );

  weights = weights ./ sum( weights(:) );
end


function [ w, nIter, flag, objValues ] = makePrecompWeights_2D_GP( traj, N, gamma, mu, varargin )

  %defaultAlg = 'fista';
  %defaultAlg = 'fista_wExtrap';
  defaultAlg = 'proxGrad_wExtrap';
  %defaultAlg = 'fista_wRestart';

  p = inputParser;
  p.addOptional( 'nIter', [], @ispositive );
  p.addParameter( 'alg', defaultAlg, @(x) true );
  p.addParameter( 'fast', true, @islogical );
  p.parse( varargin{:} );
  nIter = p.Results.nIter;
  alg = p.Results.alg;
  fast = p.Results.fast;

  if numel( nIter ) == 0
    if strcmp( alg, 'proxGrad_wExtrap' )
      nIter = 2000;
    else
      nIter = 250;
    end
  end
  if numel( alg ) == 0, alg = defaultAlg; end

  segLength = 2500;
  nTraj = size( traj, 1 );
  nSegs = ceil( nTraj / segLength );
  D = size( traj, 2 );
  if numel( N ) == 1, N = N * ones( D ); end
  if numel( gamma ) == 0, gamma = 0.25 * N; end
  if numel( gamma ) == 1, gamma = gamma * ones( numel(N), 1 ); end

  function out = applySubA( in, indxs )
    expNGammas = exp( -N ./ gamma );

    out = zeros( numel( in ), 1 );
    for t = 1 : numel( indxs )
      tIndx = indxs( t );

      nu = 2 * pi * ( bsxfun( @minus, traj, traj( tIndx, : ) ) );
      nuN = bsxfun( @times, nu, N );

      if numel( gamma ) > 0  &&  max( gamma ) < Inf
        nuGamma = bsxfun( @times, nu, gamma );

        tmp = 1 - bsxfun( @times, expNGammas, cos( nuN ) - nuGamma .* sin( nuN ) );
        factors = bsxfun( @rdivide, 2 * gamma, 1 + nuGamma.^2 ) .* tmp;
        
      else
        factors = 2 * sin( nuN ) ./ nu;
        for dIndx = 1 : D
          dFactors = factors( :, dIndx );
          dFactors( nu(:,dIndx) == 0 ) = 2 * N( dIndx );
          factors( :, dIndx ) = dFactors;
        end
      end

      colA = 2 * prod( factors, 2 );
      out = out + ( colA * in( tIndx ) );
    end
  end

  function out = gGradHat( w, subIndxs )
    if numel( subIndxs ) == 0
      out = numel( w );
      return;
    end

    out = applySubA( w, subIndxs );
  end


  function out = applyA( in, op )
    if nargin < 2, op = 'notransp'; end

    nIn = numel( in );
    out = cell( nIn, 1 );

    expNGammas = exp( -N ./ gamma );

    %p = parforProgress( nIn );
    parfor tIndx = 1 : nIn
      %p.progress( tIndx, 2500 );

      rowInSum = 0;
      eIndx = 0;
      for segIndx = 1 : nSegs
        sIndx = eIndx + 1;  % start index
        eIndx = min( segIndx * segLength, nTraj );  % end index

        prodRowIn = 2 * in( sIndx : eIndx );   %#ok<PFBNS>

        if strcmp( op, 'notransp' )  % notransp
          nu = 2 * pi * bsxfun( @minus, traj( tIndx, : ), traj( sIndx : eIndx, : ) );   %#ok<PFBNS>
        else
          nu = 2 * pi * bsxfun( @minus, traj( sIndx : eIndx, : ), traj( tIndx, : ) );
        end

        nuN = bsxfun( @times, nu, N );

        if numel( gamma ) > 0  &&  max( gamma ) < Inf
          nuGamma = bsxfun( @times, nu, gamma );

          tmp = 1 - bsxfun( @times, expNGammas, cos( nuN ) -  nuGamma .* sin( nuN ) );
          factors = bsxfun( @rdivide, 2 * gamma, 1 + nuGamma.^2 ) .* tmp;

        else
          factors = 2 * sin( nuN ) ./ nu;
          for dIndx = 1 : D
            dFactors = factors( :, dIndx );
            dFactors( nu( :, dIndx ) == 0 ) = 2 * N( dIndx );
            factors( :, dIndx ) = dFactors;
          end
        end

        prodRowIn = prodRowIn .* prod( factors, 2 );

        rowInSum = rowInSum + sum( prodRowIn );
      end

      out{ tIndx } = rowInSum;
    end
    %p.clean;

    out = cell2mat( out );
  end


  function out = makeA()

    expNGammas = exp( -N ./ gamma );

    out = zeros( nTraj, nTraj );

    parfor tIndx = 1 : nTraj
      %if mod( tIndx, 1000 ) == 0, disp( [ 'Iteration ', num2str( tIndx ) ] ); end

      tRow = zeros( 1, nTraj );
      eIndx = 0;
      for segIndx = 1 : nSegs
        sIndx = eIndx + 1;  % start index
        eIndx = min( segIndx * segLength, nTraj );  % end index

        nu = 2 * pi * bsxfun( @minus, traj( tIndx, : ), traj( sIndx : eIndx, : ) );   %#ok<PFBNS>
        nuN = bsxfun( @times, nu, N );

        if numel( gamma ) > 0  &&  max( gamma ) < Inf
          nuGamma = bsxfun( @times, nu, gamma );

          tmp = 1 - bsxfun( @times, expNGammas, cos( nuN ) -  nuGamma .* sin( nuN ) );
          factors = bsxfun( @rdivide, 2 * gamma, 1 + nuGamma.^2 ) .* tmp;

        else
          factors = 2 * sin( nuN ) ./ nu;
          for dIndx = 1 : D
            dFactors = factors( :, dIndx );
            dFactors( nu(:,dIndx) == 0 ) = 2 * N( dIndx );
            factors( :, dIndx ) = dFactors;
          end
        end

        tRow( 1, sIndx : eIndx ) = 2 * prod( factors, 2 );
        %out( tIndx, sIndx : eIndx ) = 2 * prod( factors, 2 );
      end

      out( tIndx, : ) = tRow;
    end
  end


  doAdjointCheck = false;
  if doAdjointCheck == true && ...
    ~checkAdjoint( rand( nTraj, 1 ), @applyA, 'y', rand( nTraj, 1 ) )
    error( 'A adjoint is incorrect' );
  end

  function out = g( w )
    psf = grid_2D( w, traj, 2*N, ones( size( w ) ) );
    [out,costPerPixel] = calculateCost_GP( psf, w, gamma, mu );   %#ok<ASGLU>
  end

  h = @(w) indicatorFunction( sum( w(:) ), [0.999 1.001] );

  if fast == true, A = makeA(); end

  function out = gGrad( w )
    if fast == true
      out = A * w;
    else
      out = applyA( w );
    end

    if mu > 0
      out = out + ( mu / nTraj ) * w;
    end
  end
  
  proxth = @(x,t) projectOntoProbSimplex( x );

pid = feature('getpid');
load( ['datacase_', num2str(pid), '.mat'], 'datacase' );
distWeightScaling = gamma(1) / N(1);

saveDir = [ './out/case_', indx2str(datacase,99), '/saves/scaling_', num2str(distWeightScaling) ];
algDir = [ saveDir, '/', alg ];
mkdir( saveDir );  mkdir( algDir );

  w0 = makePrecompWeights_2D( traj, 'sImg', N, 'alg', 'VORONOI' );
  %w0 = (1/size(traj,1)) * ones( size(traj,1), 1 );
  %w0 = makePrecompWeights_2D_JACKSON( traj, N );

nFile = [ saveDir, '/normA.mat' ];
if exist( nFile, 'file' )
  load( nFile, 'normA' );
else
  if fast == true
    normA = powerIteration( A, w0, 'maxIters', 10, 'verbose', true );
  else
    normA = powerIteration( @applyA, w0, 'maxIters', 10, 'verbose', true );
  end
  save( nFile, 'normA' );
end
extrapolated = [];
  grdNrm = normA + ( 0.5 * mu / nTraj );  % bound on norm of gradient
  stepSize = 0.99 / grdNrm;

  relDiffs = [];
  objValues = [];
  restarts = [];
  tol = 1d-4;
  if strcmp( 'adaptiveAcceleration', alg )
    [ w, objValues, relDiffs ] = adaptiveAccelerationOptimization( w0, @gGrad, 'proxth', proxth, ...
      'g', @g, 'h', h, 'N', nIter, 't', stepSize, 'verbose', true );

  elseif strcmp( 'fista', alg )
    [ w, objValues, relDiffs ] = fista( w0, @gGrad, proxth, 't', stepSize, 'N', nIter, ...
      'g', @g, 'h', h, 'verbose', true );

  elseif strcmp( 'fista_wExtrap', alg )
    [ w, objValues, relDiffs ] = fista_wExtrap( w0, @gGrad, proxth, 't', stepSize, 'N', nIter, ...
      'g', @g, 'h', h, 'verbose', true );

  elseif strcmp( 'fista_wLS', alg )
    t0 = 10 * stepSize;
    tol = 1d-4;
    [ w, objValues, relDiffs ] = fista_wLS( w0, @g, @gGrad, proxth, 'minStep', stepSize, ...
      'N', nIter, 't0', t0, 'tol', tol, 'h', h, 'gradNorm', grdNrm, 'verbose', true );

  elseif strcmp( 'fista_wAdaptiveRestartGradient', alg )
    [ w, objValues, relDiffs, restarts ] = fista_wAdaptiveRestartGradient( w0, @gGrad, proxth, 't', stepSize, 'N', nIter, ...
      'g', @g, 'h', h, 'tol', tol, 'verbose', true );

  elseif strcmp( 'fista_wAdaptiveRestartFunction', alg )
    [ w, objValues, relDiffs, restarts ] = fista_wAdaptiveRestartFunction( w0, @gGrad, proxth, ...
      @g, h, 't', stepSize, 'N', nIter, 'tol', tol, 'verbose', true );

  elseif strcmp( 'fista_wRestart', alg )
    q = 50;
    [ w, objValues, relDiffs, restarts ] = fista_wRestart( w0, @gGrad, proxth, q, 't', stepSize, 'N', nIter, ...
      'g', @g, 'h', h, 'tol', tol, 'verbose', true );

  elseif strcmp( 'pogm', alg )
    [ w, objValues ] = pogm( w0, @gGrad, proxth, 't', stepSize, 'N', nIter, ...
      'g', @g, 'h', h, 'verbose', true );

  elseif strcmp( 'proxGrad_wExtrap', alg )
    [ w, objValues, relDiffs, extrapolated ] = proxGrad_wExtrap( w0, @gGrad, proxth, 'N', nIter, ...
      't', stepSize, 'tol', tol, 'g', @g, 'h', h, 'verbose', true );

  elseif strcmp( 'projSubgrad', alg )
    [ w, objValues ] = projSubgrad( w0, @gGrad, @projectOntoProbSimplex, 'g', @g, ...
      't', stepSize, 'N', nIter );

  elseif strcmp( 'proxSVRG', alg )
    [ w, objValues, relDiffs ] = proxSVRG( w0, stepSize, @gGradHat, proxth, ...
      'g', @g, 'h', @h, 'nEpochs', 20, 'nStoch', 10, 'saveEvery', 10, 'verbose', true );

  elseif strcmp( 'spg', alg )
    [ w, objValues, relDiffs ] = stochasticProxGrad( w0, stepSize, @gGradHat, proxth, ...
      'g', @g, 'h', @h, 'nEpochs', 100, 'nStoch', 2000, 'saveEvery', 20, 'verbose', true );

  end
save( [ algDir, '/results.mat' ], 'w', 'objValues', 'relDiffs', 'extrapolated', 'restarts' );
  flag = 0;

  % Now find kappa scaling
  centerRegion = round( N ./ 20 );
  threshSmallK = repmat( min( 0.5 ./ N, 1d-6 ), [ size( traj, 1 ) 1 ] );
  useCenterRegion = true;
  if useCenterRegion == true
    scaledSincTraj = sin( bsxfun( @times, traj, 2 * pi * centerRegion ) ) ./ ( pi * traj );
    twoN = repmat( 2 * centerRegion, [ size( traj, 1 ) 1 ] );
  else
    scaledSincTraj = sin( bsxfun( @times, traj, 2 * pi * N ) ) ./ ( pi * traj );
    twoN = repmat( 2 * N, [ size( traj, 1 ) 1 ] );
  end
  scaledSincTraj( abs( traj ) < threshSmallK ) = twoN( abs( traj ) < threshSmallK );
  prodScaledSincTraj = prod( scaledSincTraj, 2 );
  kappa = 1 ./ ( sum( w .* prodScaledSincTraj ) );
  w = kappa * w;

  if nargout > 1, nIter = numel( objValues ); end
end


function [out,costPerPixel] = calculateCost_GP( psf, weights, distWeight, mu )
  nTraj = numel( weights );
  imgCoords = size2imgCoordinates( size( psf ) );
  [xs,ys] = meshgrid( imgCoords{2}, imgCoords{1} );
  costPerPixel = abs( psf ).^2 + 0.5 * ( mu / nTraj ) * norm( weights )^2;
  if numel( distWeight ) > 0
     costPerPixel = costPerPixel .* exp( ...
       -( ( abs(xs) / distWeight(1) ) + ( abs(ys) / distWeight(2) ) ) );
  end
  out = sum( costPerPixel(:) );
end


function weights = makePrecompWeights_2D_JACKSON( traj, N, varargin )
  % Fixed point iteration defined in "Sampling Density Compensation in MRI:
  % Rationale and an Iterative Numerical Solution" by Pipe and Menon, 1999.

  checknum = @(x) numel(x) == 0 || ( isnumeric(x) && isscalar(x) && (x >= 1) );
  p = inputParser;
  p.addParameter( 'alpha', [], checknum );
  p.addParameter( 'scaleIt', true, @islogical );
  p.addParameter( 'W', [], checknum );
  p.addParameter( 'nC', [], checknum );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  scaleIt = p.Results.scaleIt;
  W = p.Results.W;
  nC = p.Results.nC;

  % Make the Kaiser Bessel convolution kernel
  Ny = N(1);  Nx = N(2);
  Gy = Ny;
  [kCy,Cy,~] = makeKbKernel( Gy, Ny, 'alpha', alpha, 'W', W, 'nC', nC );
  Gx = Nx;
  [kCx,Cx,~] = makeKbKernel( Gx, Nx, 'alpha', alpha, 'W', W, 'nC', nC );

  nTraj = size( traj, 1 );
  weights = 1 ./ applyC_2D( ones( nTraj, 1 ), traj, traj, kCy, kCx, Cy, Cx );
  
  if scaleIt == true
    weights = weights / sum( weights(:) );
  end
end

function [weights,nIter,flag,residual] = makePrecompWeights_2D_LSQR( traj, N, varargin )
  % Optimization analog of FP

  p = inputParser;
  p.addParameter( 'alpha', [], @isnumeric );
  p.addParameter( 'W', [], @ispositive );
  p.addParameter( 'nC', [], @ispositive );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;

  Ny = 2 * N(1);   [kCy,Cy,~] = makeKbKernel( Ny, Ny, 'alpha', alpha, 'W', W, 'nC', nC );
  Nx = 2 * N(2);   [kCx,Cx,~] = makeKbKernel( Nx, Nx, 'alpha', alpha, 'W', W, 'nC', nC );

  CyCx = Cy(:) * Cx(:)';
  dkCy = zeros( size( kCy ) );  dkCy(1) = kCy(2);  dkCy(2:end) = kCy(2:end) - kCy(1:end-1);
  dkCx = zeros( size( kCx ) );  dkCx(1) = kCx(2);  dkCx(2:end) = kCx(2:end) - kCx(1:end-1);
  dkCykCx = dkCy(:) * dkCx(:)';
  scaling = 1 / sum( CyCx(:) .* dkCykCx(:) );
  Cy = Cy / scaling;
  Cx = Cx / scaling;

  nTraj = size( traj, 1 );
  C = makeC_2D( traj, traj, kCy, kCx, Cy, Cx );
  w0 = 1 ./ ( C * ones( nTraj, 1 ) );

  doCheckAdjoint = false;
  if doCheckAdjoint == true
    [checkAdjA,errCheckAdjA] = checkAdjoint( w0, @applyA );
    if checkAdjA == false
      error([ 'Check of A Adjoint failed with error: ', num2str( errCheckAdjA ) ]);
    else
      disp( 'Check of A Adjoint passed' );
    end
  end

  b = ones( size( w0 ) );
  [ weights, flag, residual, nIter ] = lsqr( C, b, [], 100, [], [], w0 );

  weights = weights / sum( weights );
end

function [weights,nIter,flag,residual] = makePrecompWeights_2D_LSQR_slow( traj, N, varargin )
  % Optimization analog of FP

  p = inputParser;
  p.addParameter( 'alpha', [], @isnumeric );
  p.addParameter( 'W', [], @ispositive );
  p.addParameter( 'nC', [], @ispositive );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;

  Ny = 2 * N(1);   [kCy,Cy,~] = makeKbKernel( Ny, Ny, 'alpha', alpha, 'W', W, 'nC', nC );
  Nx = 2 * N(2);   [kCx,Cx,~] = makeKbKernel( Nx, Nx, 'alpha', alpha, 'W', W, 'nC', nC );

  CyCx = Cy(:) * Cx(:)';
  dkCy = zeros( size( kCy ) );  dkCy(1) = kCy(2);  dkCy(2:end) = kCy(2:end) - kCy(1:end-1);
  dkCx = zeros( size( kCx ) );  dkCx(1) = kCx(2);  dkCx(2:end) = kCx(2:end) - kCx(1:end-1);
  dkCykCx = dkCy(:) * dkCx(:)';
  scaling = 1 / sum( CyCx(:) .* dkCykCx(:) );
  Cy = Cy / scaling;
  Cx = Cx / scaling;

  function out = applyA( in, op )
    if nargin < 2 || strcmp( op, 'notransp' )
      out = applyC_2D( in, traj, traj, kCy, kCx, Cy, Cx );
    elseif strcmp( op, 'transp' )
      out = applyCT_2D( in, traj, traj, kCy, kCx, Cy, Cx );
    else
      error( 'Unrecognized operator for A' );
    end
    out = out(:);
  end

  %w0 = ( 1 / size(traj,1) ) * ones( size(traj,1), 1 );
  w0 = makePrecompWeights_2D_JACKSON( traj, N, 'scaleIt', false );

  doCheckAdjoint = false;
  if doCheckAdjoint == true
    [checkAdjA,errCheckAdjA] = checkAdjoint( w0, @applyA );
    if checkAdjA == false
      error([ 'Check of A Adjoint failed with error: ', num2str( errCheckAdjA ) ]);
    else
      disp( 'Check of A Adjoint passed' );
    end
  end

  b = ones( size( w0 ) );
  [ weights, flag, residual, nIter ] = lsqr( @applyA, b, [], 100, [], [], w0 );

  weights = weights / sum( weights );
end


function [weights,lsFlag,lsRes] = makePrecompWeights_2D_SAMSANOV( ...
  kTraj, N, varargin )

  defaultAlpha = 1.5;
  checknum = @(x) numel(x) == 0 || ( isnumeric(x) && isscalar(x) && (x > 1) );
  p = inputParser;
  p.addParameter( 'alpha', defaultAlpha, checknum );
  p.addParameter( 'W', [], checknum );
  p.addParameter( 'nC', [], checknum );
  p.addParameter( 'verbose', false, @(x) islogical(x) || isnumeric(x) );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;
  verbose = p.Results.verbose;

  kbN = 2 * N;
  nGrid = ceil( alpha * kbN );
  trueAlpha = max( nGrid ./ kbN );

  % Make the Kaiser Bessel convolution kernel
  Ny=kbN(1);  Nx=kbN(2);
  Gy = Ny;
  [kCy,Cy,~] = makeKbKernel( Gy, Ny, 'alpha', trueAlpha, 'W', W, 'nC', nC );
  Gx = Nx;
  [kCx,Cx,~] = makeKbKernel( Gx, Nx, 'alpha', trueAlpha, 'W', W, 'nC', nC );

  iteration = 0;
  function out = applyA( in, type )
    if nargin > 1 && strcmp( type, 'transp' )
      reshaped = reshape( in, nGrid );
      out = applyCT_2D( reshaped, kTraj, nGrid, kCy, kCx, Cy, Cx, ...
        'type', 'noCirc' );
    else
      out = applyC_2D( in, kTraj, nGrid, kCy, kCx, Cy, Cx, ...
        'type', 'noCirc' );
      out = out(:);

      iteration = iteration + 1;
      if verbose ~= 0 && mod( iteration, 5 ) == 0
        disp(['makePrecompWeights_2D_SAMSANOV working on iteration ', ...
          num2str(iteration) ]);
      end
      %residuals(iteration) = norm( out(:) - b(:), 2 ) / norm(b(:),2);
      %if iteration>10, plot( residuals(1:iteration) ); drawnow; end
    end
  end

  nTraj = size(kTraj,1);
  b = ones( prod(nGrid), 1 );
  function y = linop_4_tfocs( in, mode )
    switch mode
      case 0
        y = [numel(b), nTraj];
      case 1
        y = applyA( in );
      case 2
        y = applyA( in, 'transp' );
    end
  end
  %varargout = linop_test( @linop_4_tfocs );

  opts = tfocs;
  opts.alg = 'N83';
  opts.maxIts = 2000;
  opts.printEvery = 1;
  %opts.tol = 1d-8;
  opts.tol = 1d-10;
  w0 = ones( nTraj, 1 );
  [weights,tfocsOptOut] = tfocs( smooth_quad, { @linop_4_tfocs, -b }, ...
    proj_Rplus, w0, opts );

  if nargout > 2
    lsFlag = tfocsOptOut.niter == opts.maxIts;
  end
  if nargout > 3
    tmp = applyA( weights );
    lsRes = norm( tmp - 1, 2 );
  end
end


function fullWeights = makePrecompWeights_2D_VORONOI( fullKTraj )

  [ kTraj, ~, ic ] = unique( fullKTraj, 'rows' );
  nTraj = size( kTraj, 1 );

  % Voronoi density compensation according to "Resampling of Data Between
  % Arbitrary Grids Using Convolution Interpolation" by Rasche et al.

  [ outerConvHullIndxs, areaOuter ] = convhull( kTraj );
  outerConvHullIndxs = outerConvHullIndxs( 1 : end - 1 );
  outerConvHullIndxs = sort( outerConvHullIndxs );
  outerTraj = kTraj(outerConvHullIndxs,:);
  innerTraj = setdiff( kTraj, outerTraj, 'rows' );

  [ ~, areaInner ] = convhull( innerTraj );
  alpha = sqrt( areaOuter / areaInner );
  newOuterTraj = alpha * outerTraj;

  augTraj = [ kTraj; newOuterTraj; ];
  augWeights = calculateVoronoiAreas( augTraj );
  weights = augWeights( 1 : nTraj );

  fullWeights = weights( ic );
  if numel( fullWeights ) ~= numel( weights )
    for trajIndx = 1 : nTraj
      fullWeights( ic == trajIndx ) = fullWeights( ic == trajIndx ) / sum( ic == trajIndx );
    end
  end

  nWeightsNotFinite = sum( ~isfinite( fullWeights ) );
  nFullTraj = size( fullKTraj, 1 );
  if nWeightsNotFinite > 0
    fullWeights( ~isfinite( fullWeights ) ) = 1 / nFullTraj;
  end
end


