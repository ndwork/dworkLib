
function [kC,C,c] = makeKbKernel( G, N, varargin )
  % This function makes the Kaiser Bessel kernel; this kernel is a
  % fundamental part of the Gridding (and inverse Gridding) algorithm.
  % [kc,C,c] = makeKbKernel( G, N,  [, 'alpha', alpha, 'W', W, 'nC', nC ] )
  %
  % Inputs:
  %   G - the number of elements of the oversampled grid
  %   N - the number of elements in the function c
  %   alpha - the oversampling factor
  %   W - the width of kernel (in units of oversampled grid points)
  %   nC - the number of elements in (half) the kernel
  % Outputs:
  %   kC - The frequency values (between 0 and 0.5) of the kernel C
  %   C - The intensity values of the kernel at the frequencies of kC
  %   c - The fourier transform of the KB kernel C

  defaultAlpha = 1.5;
  defaultW = 8;
  defaultNc = 500;
  checknum = @(x) isnumeric(x) && isscalar(x) && (x >= 1);
  p = inputParser;
  p.addOptional( 'alpha', defaultAlpha, @(x) numel(x) == 0 || checknum(x) );
  p.addOptional( 'W', defaultW, @(x) numel(x) == 0 || checknum(x) );
  p.addOptional( 'nC', defaultNc, @(x) numel(x) == 0 || checknum(x) );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;

  if numel( alpha ) == 0, alpha = defaultAlpha; end
  if numel( W ) == 0, W = defaultW; end
  if numel( nC ) == 0, nC = defaultNc; end

  kw = W/G;  % kernel width (in k-space samples)
  kC = transpose(linspace(0, 0.5*kw, nC));

  beta = pi * sqrt( W*W/(alpha*alpha) * (alpha-0.5)^2 - 0.8 );
  C = 1/kw * besseli( 0, beta * sqrt( 1 - ( 2*kC/kw ).^2 ) );

  maxC = max( C );
  C = C / maxC;

  x = size2imgCoordinates( N );
  tmp = sqrt( (pi*kw*x).^2 - beta*beta );
  c = sinc( tmp / pi );
  c = c / maxC;
end
