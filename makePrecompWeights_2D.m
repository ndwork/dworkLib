
function [weights,flag,res] = makePrecompWeights_2D( kTraj, N, varargin )
  % [weights,flag,res] = makePrecompWeights_2D( kTraj, N, ...
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
  %     CLSDC (default) - Constrainted Least Squares on trajectory points
  %     FP - specifies fixed point iteration
  %     SAMSANOV - Constrainted Least Squares on grid points
  %     VORONOI
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

  defaultAlg = 'FP';
  defaultNIter = 5;
  radImg = makeRadialImg( 2 * N );
	defaultPsfMask = radImg < abs(min(N));
  checknum = @(x) numel(x) == 0 || ( isnumeric(x) && isscalar(x) && (x >= 1) );
  p = inputParser;
  p.addParameter( 'alpha', [], checknum );
  p.addParameter( 'W', [], checknum );
  p.addParameter( 'nC', [], checknum );
  p.addParameter( 'alg', defaultAlg );
  p.addParameter( 'nIter', defaultNIter, checknum );
  p.addParameter( 'psfMask', defaultPsfMask );
  p.addParameter( 'verbose', false, @(x) islogical(x) || isnumeric(x) );
  p.parse( varargin{:} );
  alg = p.Results.alg;
  alpha = p.Results.alpha;
  W = p.Results.W;
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
      % Pipe's method
      [weights,flag,res] = makePrecompWeights_2D_FP( kTraj, N, psfMask, ...
        'alpha', alpha, 'W', W, 'nC', nC, 'nIter', nIter, 'verbose', verbose );

    case 'SAMSANOV'
      % Method from Samsanov's abstract
      [weights,flag,res] = makePrecompWeights_2D_SAMSANOV(...
        kTraj, N, 'alpha', alpha, 'W', W, 'nC', nC );

    case 'voronoi'
      % Voronoi Density compensation implementation
      weights = makePrecompWeights_2D_VORONOI( kTraj, N );

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

  scale = getScaleOfPSF( weights, traj, N, 'imgTitle', 'rtbLSDC');
  weights = scale * weights;
  close;
end


function [weights,flag,res] = makePrecompWeights_2D_FP( ...
  traj, N, psfMask, varargin )
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
  % alpha specifies the transition band of the KB filter
  kbN = 2 * N;
  Ny = kbN(1);  Nx = kbN(2);
  Gy = Ny;
  [kCy,Cy,~] = makeKbKernel( Gy, Ny, 'alpha', alpha, 'W', W, 'nC', nC );
  Gx = Nx;
  [kCx,Cx,~] = makeKbKernel( Gx, Nx, 'alpha', alpha, 'W', W, 'nC', nC );

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

  if nargout > 1
    flag = 0;
  end
  if nargout > 2
    res = norm( weights - oldWeights, 2 ) / norm( weights, 2 );
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


function [weights,flag,res] = makePrecompWeights_2D_VORONOI( kTraj, N )
  % Voronoi density compensation

  weights = voronoidens( kTraj );
	weights = min( weights, 1/size(kTraj,1) );

  radialImg = makeRadialImg( 2*N );
  mask = double( radialImg <= min(N) );

  scale = getScaleOfPSF( weights, kTraj, N, mask, 'imgTitle', 'FP' );
  weights = scale * weights;
  close
  
  flag = 0;  res = -1;
end





function [scale,mse] = getScaleOfPSF( weights, traj, N, varargin )
  % scale = getScaleOfPSF( weights, traj, N, mask [, ...
  %   'savefile', savefile, 'imgTitle', imgTitle ] )
  %
  % Written by Nicholas Dwork - Copyright 2016

  p = inputParser;
  p.addOptional( 'mask', [] );
  p.addParameter( 'savefile', [] );
  p.addParameter( 'imgTitle', [] );
  p.parse( varargin{:} );
  mask = p.Results.mask;
  savefile = p.Results.savefile;
  imgTitle = p.Results.imgTitle;

  nGrid = 2 * N;
  W = 8;
  nC = 500;

  psf = iGridT_2D( weights, traj, nGrid, 'W', W, 'nC', nC );
  psf = cropData( psf, 2*N );
  if isempty(mask), mask = ones( size(psf) ); end
  mask = cropData( mask, nGrid );
  psf = mask .* psf;
  scale = 1 ./ max(real(psf(:)));

  %psf = psf * scale;
  %absPsf=abs(psf);
  %absPsf_dB = 20*log10(absPsf);
  %absPsf_dB(~isfinite(absPsf_dB)) = 0;

  %figure;
  %imshow( imresize(absPsf_dB,2,'nearest'), [-100 0] );
  %if ~isempty(imgTitle), title(imgTitle); end
  %drawnow;

  %if ~isempty(savefile)
  %  scaled_dB = scaleImg( absPsf_dB, [-100 0], [0 1] );
  %  imwrite( scaled_dB, savefile );
  %end
  
  if nargout > 1
    b=zeros(size(psf)); b(1,1)=1; b=fftshift(b);
    mse = sum( abs( mask(:).*psf(:) - mask(:).*b(:) ).^2 ) ./ sum(mask(:));
  end
end




