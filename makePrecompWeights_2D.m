
function [weights,nOptIter,flag,res] = makePrecompWeights_2D( kTraj, varargin )
  % [weights,nOptIter,flag,res] = makePrecompWeights_2D( kTraj [, 'sImg', sImg, ...
  %   'alpha', alpha, 'W', W, 'nC', nC, 'alg', alg ] )
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
  %     GP - gradient projection algorithm
  %     GP_sparse - sparse matrix approximation of GP algorithm
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
    disp( 'Usage: [weights,nOptIter,flag,res] = makePrecompWeights_2D( kTraj, N [, ...' );
    disp( '  ''alpha'', alpha, ''W'', W, ''nC'', nC, ''alg'', alg ] )' );
    if nargout > 0, weights = []; end
    if nargout > 1, nOptIter = []; end
    if nargout > 2, flag = []; end
    if nargout > 3, res = []; end
    return
  end

  if ~isreal( kTraj ), kTraj = [ real( kTraj(:) ) imag( kTraj(:) ) ]; end

  defaultAlg = 'VORONOI';
  p = inputParser;
  p.addParameter( 'alpha', [], @isnumeric );
  p.addParameter( 'gamma', [], @isnumeric );
  p.addParameter( 'mu', 0, @isnumeric );
  p.addParameter( 'nC', [], @isnumeric );
  p.addParameter( 'nIter', [], @ispositive );
  p.addParameter( 'sImg', [], @ispositive );
  p.addParameter( 'W', [], @ispositive );
  p.addParameter( 'alg', defaultAlg, @(x) true );
  p.addParameter( 'verbose', false, @(x) islogical(x) || isnumeric(x) );
  p.parse( varargin{:} );
  alg = p.Results.alg;
  alpha = p.Results.alpha;
  gamma = p.Results.gamma;
  mu = p.Results.mu;
  nC = p.Results.nC;
  nIter = p.Results.nIter;
  sImg = p.Results.sImg;
  verbose = p.Results.verbose;
  W = p.Results.W;

  if numel( alg ) == 0, alg = defaultAlg; end
  if numel( sImg ) == 1, sImg = [ sImg sImg ]; end

  nOptIter = 1;
  flag = 0;
  res = 0;
  switch alg
    case 'CLS'
      % Constrained least squares density compensation
      [weights,nOptIter,flag,res] = makePrecompWeights_2D_CLS( kTraj, sImg );

    case 'FP'
      % Fixed point iteration
      [weights,flag,res] = makePrecompWeights_2D_FP( kTraj, sImg, 'verbose', verbose, ...
        'alpha', alpha, 'W', W, 'nC', nC );

    case 'GP'
      % Gradient Projection method
      [weights,nOptIter,flag,res] = makePrecompWeights_2D_GP( kTraj, sImg, gamma, mu, nIter );

    case 'GP_sparse'
      % sparse approximation to GP method
      [weights,nOptIter,flag,res] = makePrecompWeights_2D_GP_sparse( kTraj, sImg, gamma, mu, nIter );
      
    case 'JACKSON'
      weights = makePrecompWeights_2D_JACKSON( kTraj, sImg, 'alpha', alpha, 'W', W, 'nC', nC );

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


function scale = getScaleOfPSF( weights, traj, N )
  psf = grid_2D( weights, traj, 2*N, ones( size( weights ) ) );
  psf = cropData( psf, 0.1 * N );
  scale = 1 ./ sum( psf(:) );
end


function [ w, nIter, flag, objValues ] = makePrecompWeights_2D_GP_old( traj, N, gamma, mu, varargin )

  defaultNIter = 250;  % 300 works well
  p = inputParser;
  p.addOptional( 'nIter', defaultNIter, @ispositive );
  p.parse( varargin{:} );
  nIter = p.Results.nIter;

  if numel( nIter ) == 0, nIter = defaultNIter; end

  segLength = 25000;
  nTraj = size( traj, 1 );
  nSegs = ceil( nTraj / segLength );
  D = size( traj, 2 );
  if numel( N ) == 1, N = N * ones( D ); end
  if numel( gamma ) == 0, gamma = 0.25 * N; end
  if numel( gamma ) == 1, gamma = gamma * ones( numel(N), 1 ); end

  function out = applyA( in, op )
    if nargin < 2, op = 'notransp'; end

    nIn = numel( in );
    out = cell( nIn, 1 );

    expNGammas = exp( -N ./ gamma );

    parfor tIndx = 1 : nIn
      rowInSum = 0;
      eIndx = 0;
      for segIndx = 1 : nSegs
        sIndx = eIndx + 1;  % start index
        eIndx = min( segIndx * segLength, nTraj );  % end index

        prodRowIn = 2 * in( sIndx : eIndx );   %#ok<PFBNS>

        for dIndx = 1 : D
          Nd = N( dIndx );   %#ok<PFBNS>
          thisGamma = gamma( dIndx );   %#ok<PFBNS>
          expNdGamma = expNGammas( dIndx );   %#ok<PFBNS>

          if strcmp( op, 'notransp' )  % notransp
            nu = 2 * pi * ( traj( tIndx, dIndx ) - traj( sIndx : eIndx, dIndx ) );   %#ok<PFBNS>
          else
            nu = 2 * pi * ( traj( sIndx : eIndx, dIndx ) - traj( tIndx, dIndx ) );
          end

          if numel( thisGamma ) > 0  &&  thisGamma < Inf
            tmp = expNdGamma * ( cos( nu * Nd ) - thisGamma * nu .* sin( nu * Nd ) );
            factors = ( 2 * thisGamma ) ./ ( 1 + ( thisGamma * nu ).^2 ) .* ( 1 - tmp );

          else
            factors = 2 * sin( nu * Nd ) ./ nu;
            factors( nu == 0 ) = 2 * Nd;
          end

          prodRowIn = prodRowIn .* factors;
        end

        rowInSum = rowInSum + sum( prodRowIn );
      end

      out{ tIndx } = rowInSum;
    end
    out = cell2mat( out );
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

  h = @(w) 0;

  function out = gGrad( w )
    out = applyA( w );

    if mu > 0
      out = out + ( mu / nTraj ) * w;
    end
  end

  proxth = @(x,t) projectOntoProbSimplex( x );

  %w0 = makePrecompWeights_2D( traj, N, 'alg', 'VORONOI' );
  %w0 = (1/size(traj,1)) * ones( size(traj,1), 1 );
  w0 = makePrecompWeights_2D_JACKSON( traj, N, 'scaleIt', false );

pid = feature('getpid');
dFile = ['datacase_', num2str(pid), '.mat'];
load( dFile, 'datacase' );

nFile = ['normA_', indx2str(datacase,99), '.mat'];
if exist( nFile, 'file' )
  load( nFile, 'normA' );
else
  normA = powerIteration( @applyA, w0, 'maxIters', 10, 'verbose', true );
  save( nFile, 'normA' );
end
  grdNrm = normA + ( 0.5 * mu / nTraj );  % bound on norm of gradient
  stepSize = 0.99 / grdNrm;

  alg = 'fista_wLS';
  if strcmp( 'fista_wLS', alg )
    t0 = 10 * stepSize;
    tol = 1d-3;
    [w,objValues,relDiffs] = fista_wLS( w0, @g, @gGrad, proxth, 'minStep', stepSize, ...
      'N', nIter, 't0', t0, 'tol', tol, 'h', h, 'gradNorm', grdNrm, 'verbose', true );

save( [ 'relDiffs_', indx2str( datacase, 99 ) ], 'w', 'objValues', 'relDiffs', 'normA' );

  elseif strcmp( 'pogm', alg )
    stepSize = 0.999 / grdNrm;
    [w,objValues] = pogm( w0, @gGrad, proxth, 't', stepSize, 'N', nIter, ...
      'g', @g, 'h', h, 'verbose', true );

  elseif strcmp( 'projSubgrad', alg )
    [w,objValues] = projSubgrad( w0, @gGrad, @projectOntoProbSimplex, 'g', @g, 't', stepSize, 'N', nIter );
  end 
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


function [ w, nIter, flag, objValues ] = makePrecompWeights_2D_GP( traj, N, gamma, mu, varargin )

  defaultNIter = 250;  % 300 works well
  p = inputParser;
  p.addOptional( 'nIter', defaultNIter, @ispositive );
  p.parse( varargin{:} );
  nIter = p.Results.nIter;

  if numel( nIter ) == 0, nIter = defaultNIter; end

  segLength = 25000;
  nTraj = size( traj, 1 );
  nSegs = ceil( nTraj / segLength );
  D = size( traj, 2 );
  if numel( N ) == 1, N = N * ones( D ); end
  if numel( gamma ) == 0, gamma = 0.25 * N; end
  if numel( gamma ) == 1, gamma = gamma * ones( numel(N), 1 ); end

  function out = applyA( in, op )
    if nargin < 2, op = 'notransp'; end

    nIn = numel( in );
    out = cell( nIn, 1 );

    expNGammas = exp( -N ./ gamma );

    parfor tIndx = 1 : nIn
      rowInSum = 0;
      eIndx = 0;
      for segIndx = 1 : nSegs
        sIndx = eIndx + 1;  % start index
        eIndx = min( segIndx * segLength, nTraj );  % end index

        prodRowIn = 2 * in( sIndx : eIndx );   %#ok<PFBNS>

        if strcmp( op, 'notransp' )  % notransp
          nu = 2 * pi * ( traj( tIndx, : ) - traj( sIndx : eIndx, : ) );   %#ok<PFBNS>
        else
          nu = 2 * pi * ( traj( sIndx : eIndx, : ) - traj( tIndx, : ) );
        end

        nuNd = bsxfun( @times, nu, N );
        if numel( gamma ) > 0  &&  max( gamma ) < Inf
          tmp = bsxfun( @times, expNGammas, ...
            cos( nuNd ) - bsxfun( @times, nu .* sin( nuNd ), gamma ) );
          factors = bsxfun( @rdivide, 2 * gamma, ...
            1 + ( bsxfun( @times, nu, gamma ).^2 ) .* ( 1 - tmp ) );

        else
          factors = 2 * sin( nuNd ) ./ nu;
          for dIndx = 1 : D
            dFactors = factors( :, dIndx );
            dFactors( nu(:,dIndx) == 0 ) = 2 * N( dIndx );
            factors( :, dIndx ) = dFactors;
          end
        end

        prodRowIn = prodRowIn .* prod( factors, 2 );

        rowInSum = rowInSum + sum( prodRowIn );
      end

      out{ tIndx } = rowInSum;
    end
    out = cell2mat( out );
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

  h = @(w) 0;

  function out = gGrad( w )
    out = applyA( w );

    if mu > 0
      out = out + ( mu / nTraj ) * w;
    end
  end

  proxth = @(x,t) projectOntoProbSimplex( x );

  %w0 = makePrecompWeights_2D( traj, N, 'alg', 'VORONOI' );
  %w0 = (1/size(traj,1)) * ones( size(traj,1), 1 );  
  w0 = makePrecompWeights_2D_JACKSON( traj, N, 'scaleIt', false );

pid = feature('getpid');
dFile = ['datacase_', num2str(pid), '.mat'];
load( [ './saves/mu_', num2str(mu), '/', dFile ], 'datacase' );

nFile = [ './saves/mu_', num2str(mu), '/normA_', indx2str(datacase,99), '.mat'];
if exist( nFile, 'file' )
  load( nFile, 'normA' );
else
  normA = powerIteration( @applyA, w0, 'maxIters', 10, 'verbose', true );
  save( nFile, 'normA' );
end
  grdNrm = normA + ( 0.5 * mu / nTraj );  % bound on norm of gradient
  stepSize = 0.99 / grdNrm;

  alg = 'fista_wLS';
  if strcmp( 'fista_wLS', alg )
    t0 = 10 * stepSize;
    tol = 1d-3;
    [w,objValues,relDiffs] = fista_wLS( w0, @g, @gGrad, proxth, 'minStep', stepSize, ...
      'N', nIter, 't0', t0, 'tol', tol, 'h', h, 'gradNorm', grdNrm, 'verbose', true );

save( [ './saves/mu_', num2str(mu), '/', 'relDiffs_', indx2str( datacase, 99 ) ], ...
  'w', 'objValues', 'relDiffs', 'normA' );

  elseif strcmp( 'pogm', alg )
    stepSize = 0.999 / grdNrm;
    [w,objValues] = pogm( w0, @gGrad, proxth, 't', stepSize, 'N', nIter, ...
      'g', @g, 'h', h, 'verbose', true );

  elseif strcmp( 'projSubgrad', alg )
    [w,objValues] = projSubgrad( w0, @gGrad, @projectOntoProbSimplex, 'g', @g, 't', stepSize, 'N', nIter );
  end 
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



function [ w, nIter, flag, objValues ] = makePrecompWeights_2D_GP_sparse( traj, N, gamma, mu, varargin )
  
  defaultNIter = 100;
  p = inputParser;
  p.addOptional( 'nIter', defaultNIter, @ispositive );
  p.parse( varargin{:} );
  nIter = p.Results.nIter;

  if numel( nIter ) == 0, nIter = defaultNIter; end

  nTraj = size( traj, 1 );
  D = size( traj, 2 );
  if numel( N ) == 1, N = N * ones( D ); end
  if numel( gamma ) == 0, gamma = 0.25 * N; end
  if numel( gamma ) == 1, gamma = gamma * ones( numel(N), 1 ); end


pid = feature('getpid');
dFile = ['datacase_', num2str(pid), '.mat'];
load( dFile, 'datacase' );


aSparseFile = ['A_sparse_', indx2str(datacase,99), '.mat'];

if exist( aSparseFile, 'file' )
  load( aSparseFile, 'A' );
else
  A = makeSparseA( traj, N, gamma, 2.0 );
  save( aSparseFile, 'A' );
end


  function out = applyA( in, op )
    if nargin < 2 || strcmp( op, 'notransp' )
      out = A * in;
    else
      out = A' * in;
    end
  end

  function out = g( w )
    psf = grid_2D( w, traj, 2*N, ones( size( w ) ) );
    [out,costPerPixel] = calculateCost_GP( psf, w, gamma, mu );   %#ok<ASGLU>
  end

  h = @(w) 0;

  function out = gGrad( w )
    out = A * w;

    if mu > 0
      out = out + ( mu / nTraj ) * w;
    end
  end

  proxth = @(x,t) projectOntoProbSimplex( x );



  %w0 = makePrecompWeights_2D( traj, N, 'alg', 'VORONOI' );
  w0 = (1/size(traj,1)) * ones( size(traj,1), 1 );
  %w0 = makePrecompWeights_2D_JACKSON( traj, N, 'scaleIt', false );


nFile = ['normA_sparse_', indx2str(datacase,99), '.mat'];
if exist( nFile, 'file' )
  load( nFile, 'normA' );
else
  normA = powerIteration( @applyA, w0, 'maxIters', 10, 'verbose', true );
  save( nFile, 'normA' );
end
  grdNrm = normA + ( 0.5 * mu / nTraj );  % bound on norm of gradient
  stepSize = 0.99 / grdNrm;

  alg = 'fista_wLS';
  if strcmp( 'fista_wLS', alg )
    t0 = 10 * stepSize;
    tol = 1d-3;
    [w,objValues,relDiffs] = fista_wLS( w0, @g, @gGrad, proxth, 'minStep', stepSize, ...
      'N', nIter, 't0', t0, 'tol', tol, 'h', h, 'gradNorm', grdNrm, 'verbose', true );

%load( 'datacase.mat', 'datacase' );
save( [ 'relDiffs_sparse_', indx2str( datacase, 99 ) ], 'w', 'objValues', 'relDiffs', 'normA' );
    
  elseif strcmp( 'pogm', alg )
    stepSize = 0.999 / grdNrm;
    [w,objValues] = pogm( w0, @gGrad, proxth, 't', stepSize, 'N', nIter, ...
      'g', @g, 'h', h, 'verbose', true );

  elseif strcmp( 'projSubgrad', alg )
    [w,objValues] = projSubgrad( w0, @gGrad, @projectOntoProbSimplex, 'g', @g, 't', stepSize, 'N', nIter );
  end 
  %figure;  plotnice( objValues );
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


function A = makeSparseA( traj, N, gamma, thresh )

  nTraj = size( traj, 1 );
  rowIndxs = cell( nTraj, 1 );
  colIndxs = cell( nTraj, 1 );
  valuesA = cell( nTraj, 1 );
  D = numel( N );

  segLength = 30000;

  p = parforProgress( nTraj );
  parfor trajIndx = 1 : nTraj
    p.progress( trajIndx, 250 );   %#ok<PFBNS>
    a = 2 * ones( nTraj, 1 );

    for dIndx = 1 : D
      Nd = N( dIndx );   %#ok<PFBNS>
      thisGamma = gamma( dIndx );   %#ok<PFBNS>

      eIndx = 0;
      nSegs = ceil( nTraj / segLength );
      for segIndx = 1 : nSegs
        sIndx = eIndx + 1;  % start index
        eIndx = min( segIndx * segLength, nTraj );  % end index

        nu = 2 * pi * ( traj( trajIndx, dIndx ) - traj( sIndx:eIndx, dIndx ) );   %#ok<PFBNS>

        if numel( thisGamma ) > 0  &&  thisGamma < Inf
          tmp1 = ( -thisGamma * exp( -Nd / thisGamma ) ) * cos( nu * Nd ) ;
          tmp2 = ( thisGamma * thisGamma * exp( -Nd / thisGamma ) ) * nu .* sin( nu * Nd );
          tmp = tmp1 + tmp2 + thisGamma;
          factors = 2 ./ ( 1 + ( thisGamma * nu ).^2 ) .* tmp;
          factors( nu == 0 ) = 2 * thisGamma * ( 1 - exp( -Nd / thisGamma ) );

        else
          factors = 2 * sin( nu * Nd ) ./ nu;
          factors( nu == 0 ) = 2 * Nd;
        end

        a( sIndx : eIndx ) = a( sIndx : eIndx ) .* factors;
      end
    end

    rowIndxs{ trajIndx } = trajIndx * ones( sum( a > thresh ), 1 );
    colIndxs{ trajIndx } = find( a > thresh );
    valuesA{ trajIndx } = a( a > thresh );
  end
  p.clean;

  rowIndxs = cell2mat( rowIndxs );
  colIndxs = cell2mat( colIndxs );
  valuesA = cell2mat( valuesA );

  A = sparse( colIndxs, rowIndxs, valuesA, nTraj, nTraj );
end





function weights = makePrecompWeights_2D_JACKSON( traj, N, varargin )
  % Fixed point iteration defined in "Sampling Density Compensation in MRI:
  % Rationale and an Iterative Numerical Solution" by Pipe and Menon, 1999.

  checknum = @(x) numel(x) == 0 || ( isnumeric(x) && isscalar(x) && (x >= 1) );
  p = inputParser;
  p.addParameter( 'alpha', [], checknum );
  p.addParameter( 'W', [], checknum );
  p.addParameter( 'nC', [], checknum );
  p.addParameter( 'scaleIt', true, @islogical );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;
  scaleIt = p.Results.scaleIt;

  % Make the Kaiser Bessel convolution kernel
  Ny = N(1);  Nx = N(2);
  Gy = Ny;
  [kCy,Cy,~] = makeKbKernel( Gy, Ny, 'alpha', alpha, 'W', W, 'nC', nC );
  Gx = Nx;
  [kCx,Cx,~] = makeKbKernel( Gx, Nx, 'alpha', alpha, 'W', W, 'nC', nC );

  nTraj = size( traj, 1 );
  weights = 1 ./ applyC_2D( ones( nTraj, 1 ), traj, traj, kCy, kCx, Cy, Cx );

  if scaleIt == true
    weights = weights / sum( weights );
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
  [ weights, flag, residual, nIter ] = lsqr( @applyA, b, [], 60, [], [], w0 );

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

  [ outerConvHullIndxs, areaOuter ] = convhull( kTraj, 'Simplify', true );
  outerConvHullIndxs = outerConvHullIndxs( 1 : end - 1 );
  outerConvHullIndxs = sort( outerConvHullIndxs );
  outerTraj = kTraj(outerConvHullIndxs,:);
  innerTraj = setdiff( kTraj, outerTraj, 'rows' );

  [ ~, areaInner ] = convhull( innerTraj, 'Simplify', true );
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

end


