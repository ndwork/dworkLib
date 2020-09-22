
function out = iGridT_2D( F, traj, N, varargin )
  % out = iGridT_2D( F, traj, N, [ 'alpha', alpha, 'W', W, 'nC', nC ] )
  %
  % Gridding (without density correction) is the adjoint of MRI encoding
  % with inverse gridding.  This function applies the transpose of
  % inverse gridding to the input data.
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
  % Written by Nicholas Dwork (c) 2016
  % Based on EE369C notes written by John Pauly

  if numel(N)==1
    Ny=N;  Nx=N; 
  else
    Ny=N(1); Nx=N(2);
  end

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

  % Make the Kaiser Bessel convolution kernel
  Gy = Ny;
  [kCy,Cy,cImgY] = makeKbKernel( Gy, Ny, 'alpha', alpha, 'W', W, 'nC', nC );
  Gx = Nx;
  [kCx,Cx,cImgX] = makeKbKernel( Gx, Nx, 'alpha', alpha, 'W', W, 'nC', nC );

  % Perform a circular convolution
  fftGridded = applyC_2D( F, traj, [Ny Nx], kCy, kCx, Cy, Cx );

  % Perform an inverse fft
  data = fftshift( uifft2( ifftshift(fftGridded) ) );

  % Perform deapodization
  cImg = cImgY * transpose(cImgX);
  out = data ./ cImg;
end

