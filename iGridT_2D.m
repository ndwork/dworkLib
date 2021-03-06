
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
  %   traj is a Mx2 array specifying the k-space trajectory.
  %     The first/second column is kx/ky
  %     The units are normalized to [-0.5,0.5).
  %   N is the size of the output image
  %     If N is a scalar, then the final image is assumed to be square
  %     If N is a 2 element array, then N = [Ny Nx]
  %
  % Optional Inputs:
  %   alpha - a float parameter specifying the oversampling factor
  %   W - an integer specifying the kernel's width
  %   nC - specifies the number of samples in the kernel
  %   verbose - if set to true, outputs processing status
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.
  

  if numel(N)==1
    Ny = N;     Nx = N; 
  else
    Ny = N(1);  Nx = N(2);
  end

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
  Gy = Ny;
  [kCy,Cy,cImgY] = makeKbKernel( Gy, Ny, 'alpha', alpha, 'W', W, 'nC', nC );
  Gx = Nx;
  [kCx,Cx,cImgX] = makeKbKernel( Gx, Nx, 'alpha', alpha, 'W', W, 'nC', nC );

  % Perform a circular convolution
  fftGridded = applyC_2D( F, traj, [Ny Nx], kCy, kCx, Cy, Cx );

  % Perform an inverse fft
  data = fftshift( fftshift( ifft( ifft( ifftshift( ifftshift( ...
    fftGridded, 1 ), 2 ), [], 1 ), [], 2 ), 1 ), 2 );

  % Perform deapodization
  deapod = ( Nx * Ny ).^2 ./ ( cImgY * transpose( cImgX ) ) ;
  out = bsxfun( @times, data, deapod );
end

