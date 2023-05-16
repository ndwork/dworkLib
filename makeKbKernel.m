
function [kC,C,c] = makeKbKernel( G, varargin )
  % [kc,C,c] = makeKbKernel( G, [, N,  'alpha', alpha, 'W', W, 'nC', nC ] )
  %
  % This function makes the Kaiser Bessel kernel; this kernel is a
  % fundamental part of the Gridding (and inverse Gridding) algorithm.
  % Consider reviewing "Selection of a Convolution Function for Fourier
  % Inversion Using Gridding" by Jackson et al. for details.
  %
  % Inputs:
  %   G - the number of elements of the oversampled grid
  %   N - the number of elements in the function c
  %
  % Optional Inputs:
  %   alpha - the oversampling factor
  %   W - the width of kernel (in units of oversampled grid points)
  %   nC - the number of elements in (half) the kernel
  %
  % Outputs:
  %   kC - The frequency values (between 0 and 0.5) of the kernel C
  %   C - The intensity values of the kernel at the frequencies of kC
  %   c - The fourier transform of the KB kernel C
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  defaultAlpha = 1.5;
  defaultW = 8;
  defaultNc = 500;
  
  p = inputParser;
  checknum = @(x) numel(x) == 0 || ( isnumeric(x) && isscalar(x) && (x >= 1) );
  p.addOptional( 'N', [], checknum );
  p.addParameter( 'alpha', 1.5, checknum );
  p.addParameter( 'W', 8, checknum );
  p.addParameter( 'nC', 500, checknum );
  p.parse( varargin{:} );
  N = p.Results.N;
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;

  if numel( N ) == 0, N = G; end
  if numel( alpha ) == 0, alpha = defaultAlpha; end
  if numel( W ) == 0, W = defaultW; end
  if numel( nC ) == 0, nC = defaultNc; end

  kw = W / G;  % kernel width (in k-space samples)
  kwInv = G / W;
  kC = linspace( 0, 0.5 * kw, nC );
  kC = kC(:);

  beta = pi * sqrt( W*W / (alpha*alpha) * (alpha-0.5)^2 - 0.8 );
  C = kwInv * besseli( 0, beta * sqrt( 1 - ( 2 * kC * kwInv ).^2 ) );

  maxC = max( C(:) );
  C = C ./ maxC;

  if nargout > 2
    x = size2imgCoordinates( ceil( alpha * N ) );
    tmp = sqrt( ( pi * pi * kw * kw ) * ( x .* x ) - ( beta * beta ) );
    c = sinc( tmp / pi );
    c = c ./ maxC;
  end
end
