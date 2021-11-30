
function out = iGridT_2D( F, traj, N, varargin )
  % out = iGridT_2D( F, traj, N, [ 'alpha', alpha, 'W', W, 'nC', nC ] )
  %
  % Gridding (without density correction) is the adjoint of MRI encoding
  % with inverse gridding.  This function applies the transpose of
  % inverse gridding to the input data.
  % Based on EE369C notes written by John Pauly
  % Detailed in the following document http://nicholasdwork.com/tutorials/dworkGridding.pdf
  %
  % Inputs:
  %   F is a 1D array of M elements specifying the k-space data values
  %   traj is a Mx2 array specifying the k-space trajectory or an M element complex array.
  %     If Mx2, then the first/second column is kx/ky.  Otherwise, kx/ky is real/imag.
  %     The units are normalized to [-0.5,0.5).
  %   N is the size of the output image
  %     If N is a scalar, then the final image is assumed to be square
  %     If N is a 2 element array, then N = [Ny Nx]
  %
  % Optional Inputs:
  %   alpha - a float parameter specifying the oversampling factor
  %   W - an integer specifying the kernel's width
  %   nC - specifies the number of samples in the kernel
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
    disp( 'Usage:  out = iGridT_2D( F, traj, N, [ ''alpha'', alpha, ''W'', W, ''nC'', nC ] )' );
    if nargout > 0, out = []; end
    return;
  end

  if numel(N) == 1
    Ny = N;     Nx = N; 
  else
    Ny = N(1);  Nx = N(2);
  end

  checknum = @(x) numel(x) == 0 || min( isnumeric(x) & isscalar(x) & (x >= 1) ) == 1;
  p = inputParser;
  p.addParameter( 'alpha', [], checknum );
  p.addParameter( 'W', [], checknum );
  p.addParameter( 'nC', [], checknum );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;

  if numel( alpha ) == 0, alphaY = 1.5;  alphaX = 1.5;
  elseif numel( alpha ) == 1, alphaY = alpha;  alphaX = alpha;
  elseif numel( alpha ) == 2, alphaY = alpha(1); alphaX = alpha(2);
  end
  
  if ~isreal( traj ), traj = [ real( traj(:) )  imag( traj(:) )  ]; end

  % Make the Kaiser Bessel convolution kernel
  Gy = Ny;
  [kCy,Cy,cY] = makeKbKernel( Gy, Ny, 'alpha', alphaY, 'W', W, 'nC', nC );
  Gx = Nx;
  [kCx,Cx,cX] = makeKbKernel( Gx, Nx, 'alpha', alphaX, 'W', W, 'nC', nC );
  Nc = [ numel( cY )  numel( cX ) ];

  % Perform a circular convolution
  fftGridded = applyCT_2D( F, Nc, traj, kCy, kCx, Cy, Cx );

  % Perform an ifft
  data = fftshift2( ifft2( ifftshift2( fftGridded ) ) );

  % Perform deapodization
  out = bsxfun( @rdivide, data, cY * transpose( cX ) );

  % Crop out center region if oversampling was used
  if ( alphaY ~= 1 ) || ( alphaX ~= 1 )
    sOut = size( out );
    sOut( 1 : 2 ) = [ Ny Nx ];
    out = cropData( out, sOut );
  end
end
