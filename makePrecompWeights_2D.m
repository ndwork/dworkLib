
function [weights,flag,res] = makePrecompWeights_2D( kTraj, varargin )
  % [weights,flag,res] = makePrecompWeights_2D( kTraj [, N, ...
  %   [ 'alpha', alpha, 'W', W, 'nC', nC, 'alg', alg, 'psfMask', psfMask ] )
  %
  % Determine the density pre-compensation weights to be used in gridding
  %
  % Inputs:
  %   kTraj is a Mx2 element array specifying the k-space trajectory.
  %     The first/second column is the ky/kx location.
  %     The units are normalized to [-0.5,0.5)
  %   N is a 2 element array [Ny Nx] representing the number of grid points
  %
  % Optional Inputs:
  %   alpha - the oversampling factor > 1
  %   W - the window width in pixels
  %   nC - the number of points to sample the convolution kernel
  %   alg - a string specifying the algorithm to use
  %     FP - specifies fixed point iteration
  %     GP - gradient projection algorithm
  %     SAMSANOV - Constrainted Least Squares on grid points
  %     VORONOI (default) - uses the area of each voronoi cell as the metric of density
  %   nIter - specifies the number of iterations of fp method
  %   psfMask - only used by space domain optimizations
  %   verbose - true/false
  %
  % Outputs:
  %   weights - 1D array with density compensation weights
  %
  % Optional Outputs:
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

  defaultAlg = 'VORONOI';
  defaultNIter = 5;
  checknum = @(x) numel(x) == 0 || ( isnumeric(x) && isscalar(x) && (x >= 1) );
  p = inputParser;
  p.addOptional( 'N', [], @ispositive );
  p.addParameter( 'alpha', [], checknum );
  p.addParameter( 'gamma', 21, @isnumeric );
  p.addParameter( 'W', [], checknum );
  p.addParameter( 'mu', 0, @isnumeric );
  p.addParameter( 'nC', [], checknum );
  p.addParameter( 'alg', defaultAlg, @(x) true );
  p.addParameter( 'nIter', defaultNIter, checknum );
  p.addParameter( 'psfMask', [] );
  p.addParameter( 'verbose', false, @(x) islogical(x) || isnumeric(x) );
  p.parse( varargin{:} );
  N = p.Results.N;
  alg = p.Results.alg;
  alpha = p.Results.alpha;
  gamma = p.Results.gamma;
  W = p.Results.W;
  mu = p.Results.mu;
  nC = p.Results.nC;
  nIter = p.Results.nIter;
  psfMask = p.Results.psfMask;
  verbose = p.Results.verbose;

  flag = 0;
  res = 0;
  switch alg
    case 'CLSDC'
      % tbLSDC with non-negativity constraint
      [weights,flag,res] = makePrecompWeights_2D_CLSDC( ...
        kTraj, N, 'alpha', alpha, 'W', W, 'nC', nC );

    case 'FP'
      % Fixed point iteration
      [weights,flag,res] = makePrecompWeights_2D_FP( kTraj, N, ...
        'alpha', alpha, 'W', W, 'nC', nC, 'nIter', nIter, 'verbose', verbose );

    case 'GP'
      % Gradient Projection method
      [weights,flag,res] = makePrecompWeights_2D_GP( kTraj, N, gamma, mu );

    case 'JACKSON'
      [weights,flag] = makePrecompWeights_2D_JACKSON( kTraj, N, ...
        'alpha', alpha, 'W', W, 'nC', nC );

    case 'SAMSANOV'
      [weights,flag,res] = makePrecompWeights_2D_SAMSANOV(...
        kTraj, N, 'alpha', alpha, 'W', W, 'nC', nC );

    case 'VORONOI'
      % Voronoi Density compensation implementation
      weights = makePrecompWeights_2D_VORONOI( kTraj );

    otherwise
      error('makePrecompWeights: Algorithm not recognized');
  end

end


function [weights,flag,residual] = makePrecompWeights_2D_CLSDC( ...
  traj, N, varargin )
  % Optimization analog of FP with a non-negativity constraint

  defaultAlpha = 1.5;
  defaultW = 8;
  defaultNc = 500;
  checknum = @(x) numel(x)==0 || ( isnumeric(x) && isscalar(x) && (x > 1) );
  p = inputParser;
  p.addParameter( 'alpha', defaultAlpha, checknum );
  p.addParameter( 'W', defaultW, checknum );
  p.addParameter( 'nC', defaultNc, checknum );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;

  % Make the Kaiser Bessel convolution kernel
  kbN = 2*N;
  Ny=kbN(1);  Nx=kbN(2);
  Gy = Ny;
  [kCy,Cy,~] = makeKbKernel( Gy, Ny, 'alpha', alpha, 'W', W, 'nC', nC );
  Gx = Nx;
  [kCx,Cx,~] = makeKbKernel( Gx, Nx, 'alpha', alpha, 'W', W, 'nC', nC );

  iteration = 0;
  function out = applyA( in, type )
    if nargin > 1 && strcmp( type, 'transp' )
      out = applyCT_2D( in, traj, kCy, kCx, Cy, Cx, traj );
    else
      out = applyC_2D( in, traj, 0, kCy, kCx, Cy, Cx, traj );

      iteration = iteration + 1;
      if verbose ~= 0 && mod( iteration, 5 ) == 0
        disp(['makePrecompWeights_2D_CLSDC working on iteration ', ...
          num2str(iteration) ]);
      end
      %residuals(iteration) = norm( out(:) - b(:), 2 ) / norm(b(:),2);
      %if iteration>10, plot( residuals(1:iteration) ); drawnow; end
    end
  end

  nTraj = size(traj,1);
  b = ones(nTraj,1);
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
  opts = tfocs;
  opts.maxIts = 2000;
  opts.printEvery = 1;
  opts.tol = 1d-5;
  x0 = ones( nTraj, 1 );
  [weights,tfocsOptOut] = tfocs( smooth_quad, { @linop_4_tfocs, -b }, ...
    proj_Rplus, x0, opts );
  %weights = tfocs( smooth_huber(0.02), { @linop_4_tfocs, -b }, [], x0, opts );

  if nargout > 1
    flag = 0;  % tfocs doesn't provide flag
  end
  if nargout > 2
    psf = applyA( weights );
    residual = norm( psf(:) - b(:), 2 ) / norm(b(:),2);
  end

  close;
end


function [weights,flag,res] = makePrecompWeights_2D_FP( ...
  traj, N, varargin )
  % Fixed point iteration defined in "Sampling Density Compensation in MRI:
  % Rationale and an Iterative Numerical Solution" by Pipe and Menon, 1999.

  checknum = @(x) numel(x) == 0 || ( isnumeric(x) && isscalar(x) && (x >= 1) );
  p = inputParser;
  p.addParameter( 'alpha', [], checknum );
  p.addParameter( 'W', [], checknum );
  p.addParameter( 'nC', [], checknum );
  p.addParameter( 'nIter', 8, checknum );
  p.addParameter( 'verbose', false, @(x) islogical(x) || isnumeric(x) );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;
  nIter = p.Results.nIter;
  verbose = p.Results.verbose;

  % Make the Kaiser Bessel convolution kernel
  Ny = N(1);  Nx = N(2);
  Gy = Ny;
  [kCy,Cy,~] = makeKbKernel( Gy, Ny, 'alpha', alpha, 'W', W, 'nC', nC );
  Gx = Nx;
  [kCx,Cx,~] = makeKbKernel( Gx, Nx, 'alpha', alpha, 'W', W, 'nC', nC );

  dky = kCy(2) - kCy(1);  Cy = Cy / ( sum( Cy(:) ) / dky );
  dkx = kCx(2) - kCx(1);  Cx = Cx / ( sum( Cx(:) ) / dkx );

  nTraj = size( traj, 1 );
  weights = ones( nTraj, 1 );

  flag = 1;
  for iteration=1:nIter
    if verbose ~= 0
      disp(['makePrecompWeights_2D_FP working on iteration ', num2str(iteration) ]);
    end

    oldWeights = weights;
    denom = applyC_2D( oldWeights, traj, traj, kCy, kCx, Cy, Cx );
    weights = oldWeights ./ denom;
  end

  scale = getScaleOfPSF( weights, traj, round( 0.025 * N ) );
  weights = weights .* scale;

  %weights = weights ./ sum( weights(:) );
  
  if nargout > 1
    flag = 0;
  end
  if nargout > 2
    res = norm( weights - oldWeights, 2 ) / norm( weights, 2 );
  end
end


function scale = getScaleOfPSF( weights, traj, N )
  psf = grid_2D( weights, traj, 2*N, ones( size( weights ) ) );
  scale = 1 ./ sum( psf(:) );
end


function [ w, flag, objValues ] = makePrecompWeights_2D_GP( traj, N, gamma, mu )

  if numel( gamma ) == 1, gamma = ones( numel(N), 1 ); end

  nTraj = size( traj, 1 );
  D = size( traj, 2 );
  segLength = 30000;

  function out = applyA( in, op )
    if nargin < 2, op = 'notransp'; end

    nIn = numel( in );
    out = cell( nIn, 1 );

    parfor tIndx = 1 : nIn
      row = 2 * ones( nTraj, 1 );

      for dIndx = 1 : D
        Nd = N( dIndx );   %#ok<PFBNS>
        thisGamma = gamma( dIndx );   %#ok<PFBNS>

        eIndx = 0;
        nSegs = ceil( nTraj / segLength );
        for segIndx = 1 : nSegs
          sIndx = eIndx + 1;  % start index
          eIndx = min( segIndx * segLength, nTraj );  % end index

          if strcmp( op, 'notransp' )  % notransp
            diffKs = traj( tIndx, dIndx ) - traj( sIndx:eIndx, dIndx );   %#ok<PFBNS>
          else
            diffKs = traj( sIndx:eIndx, dIndx ) - traj( tIndx, dIndx );
          end
          nu = 2 * pi * diffKs;

          if numel( thisGamma ) > 0  &&  thisGamma < Inf
            tmp1 = ( -thisGamma * exp( -Nd / thisGamma ) ) * cos( nu * Nd ) ;
            tmp2 = ( thisGamma * thisGamma * exp( -Nd / thisGamma ) ) * nu .* sin( nu * Nd );
            tmp3 = thisGamma;
            tmp = tmp1 + tmp2 + tmp3;
            factors = 2 ./ ( 1 + ( thisGamma * nu ).^2 ) .* tmp;
            factors( nu == 0 ) = 2 * thisGamma * ( 1 - exp( -Nd / thisGamma ) );

          else
            factors = sin( nu * Nd ) ./ ( pi * diffKs );
            factors( nu == 0 ) = 2 * Nd;
          end

          row( sIndx : eIndx ) = row( sIndx : eIndx ) .* factors;
        end
      end

      out{ tIndx } = sum( row .* in );
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
  w0 = (1/size(traj,1)) * ones( size(traj,1), 1 );

  normA = powerIteration( @applyA, w0, 'maxIters', 100, 'verbose', true );
  grdNrm = normA + ( 0.5 * mu / nTraj );  % Lower bound on norm of gradient
  stepSize = 0.99 / grdNrm;

  alg = 'fista_wLS';
  if strcmp( 'fista_wLS', alg )
    stepSize = 1d-4;
    [w,objValues] = fista_wLS( w0, @g, @gGrad, proxth, ...
      'N', 60, 't0', stepSize, 'h', h, 'gradNorm', grdNrm, 'verbose', true );

  elseif strcmp( 'pogm', alg )
    [w,objValues] = pogm( w0, @gGrad, proxth, 't', stepSize, 'N', 60, ...
      'g', @g, 'h', h, 'verbose', true );

  elseif strcmp( 'projSubgrad', alg )
    [w,objValues] = projSubgrad( w0, @gGrad, @projectOntoProbSimplex, 'g', @g, 't', stepSize, 'N', 10000 );
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



function [ weights, flag ] = makePrecompWeights_2D_JACKSON( traj, N, varargin )
  % Fixed point iteration defined in "Sampling Density Compensation in MRI:
  % Rationale and an Iterative Numerical Solution" by Pipe and Menon, 1999.

  checknum = @(x) numel(x) == 0 || ( isnumeric(x) && isscalar(x) && (x >= 1) );
  p = inputParser;
  p.addParameter( 'alpha', [], checknum );
  p.addParameter( 'W', [], checknum );
  p.addParameter( 'nC', [], checknum );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
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

  weights = weights / sum( weights(:) );
  
  if nargout > 1
    flag = 0;
  end
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
  if nargout > 1
    lsFlag = tfocsOptOut.niter == opts.maxIts;
  end
  if nargout > 2
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


