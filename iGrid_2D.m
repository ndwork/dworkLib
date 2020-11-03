
function F = iGrid_2D( data, traj, varargin )
  % F = iGrid_2D( data, traj, [ 'alpha', alpha, 'W', W, 'nC', nC ] )
  %
  % Inverse Gridding based on EE369C notes by John Pauly and Beatty et. al., IEEE TMI, 2005
  % Definitions and details according to http://nicholasdwork.com/tutorials/dworkGridding.pdf
  %
  % Inputs
  %   data is a 2D array specifying the volume to be encoded
  %   traj is a Mx2 array specifying the k-space trajectory.
  %     The first/second column is ky/kx
  %     The units are normalized to [-0.5,0.5).
  %
  % Optional Inputs:
  %   alpha is the oversampling factor >= 1  (default of 1.5)
  %   W is the window width in pixels  (default of 8)
  %   nC is the number of points to sample the convolution kernel (default of 500)
  %
  % Output:
  %   F the estimates of the Fourier coefficients along the trajectory
  %
  % Written by Nicholas Dwork (c) 2016

  defaultAlpha = 1.5;
  defaultW = 8;
  defaultNc = 500;
  checknum = @(x) isnumeric(x) && isscalar(x) && (x >= 1);
  p = inputParser;
  p.addParameter( 'alpha', defaultAlpha, @(x) numel(x) == 0 || checknum(x) );
  p.addParameter( 'W', defaultW, @(x) numel(x) == 0 || checknum(x) );
  p.addParameter( 'nC', defaultNc, @(x) numel(x) == 0 || checknum(x) );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;

  if numel( alpha ) == 0, alpha = defaultAlpha; end
  if numel( W ) == 0, W = defaultW; end
  if numel( nC ) == 0, nC = defaultNc; end
  
  [Ny,Nx] = size( data );

  % Make the Kaiser Bessel convolution kernel
  Gy = Ny;
  [kCy,Cy,cImgY] = makeKbKernel( Gy, Ny, alpha, W, nC );
  Gx = Nx;
  [kCx,Cx,cImgX] = makeKbKernel( Gx, Nx, alpha, W, nC );

  % Pre-emphasize the image
  denom = cImgY * transpose(cImgX);
  preEmphasized = ( Nx * Ny ) * ( data ./ denom );

  % Perform an fft
  fftData = fftshift( fft2( ifftshift( preEmphasized ) ) );

  % Perform a circular convolution
  N = [Ny Nx];
  F = applyCT_2D( fftData, traj, N, kCy, kCx, Cy, Cx );
end

