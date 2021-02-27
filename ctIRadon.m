
function out = ctIRadon( sinogram, thetas, dSize, cx, cy, Nx, Ny, dx, dy, varargin )
  % out = ctIRadon( sinogram, thetas, dSize, cx, cy, Nx, Ny, dx, dy [, ...
  %   'window', 'Hanning/none', 'type', 'iso/fast' ] )
  %
  % Performs discrete inverse Radon transform with physical units of a Computed
  % Tomography detector
  %
  % Inputs:
  % sinogram - each row corresponds to a theta, each column to a detector
  % thetas - a 1D array of thetas corresponding to each sinogram row
  % dSize - the size of each detector in meters
  % cx - center x position of reconstructed region
  % cy - center y position of reconstructed region
  % Nx - number of pixels horizontally in reconstruction
  % Ny - number of pixels vertically in reconstruction
  % dx - the horizontal size of each pixel in meters
  % dy - the vertical size of each pixel in meters
  % window - the window to be applied to the ramp filter
  %
  % Written by Nicholas Dwork - Copyright 2013
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  defaultWin = 'none';
  expectedWins = { 'Hanning', 'none' };
  p = inputParser;
  p.addRequired( 'sinogram', @(x) ismatrix(x) );
  p.addRequired( 'thetas', @(x) ismatrix(x) );
  p.addRequired( 'dSize',@isnumeric);
  p.addRequired( 'cx', @isnumeric );
  p.addRequired( 'cy', @isnumeric );
  p.addRequired( 'Nx', @isnumeric );
  p.addRequired( 'Ny', @isnumeric );
  p.addRequired( 'dx', @isnumeric );
  p.addRequired( 'dy', @isnumeric );
  p.addOptional( 'window', defaultWin, @(x) any( validatestring(x,expectedWins) ) );
  p.parse( sinogram, thetas, dSize, cx, cy, Nx, Ny, dx, dy, varargin{:} );
  window = p.Results.window;

  [nThetas, nDetectors] = size(sinogram);

  %% Filter the sinogram
  cn = floor( 0.5 * nDetectors );
  n = ( (0:nDetectors-1) - cn );
  %dLocs = n * dSize;

  h = zeros( 1, numel(n) );
  oddNs = find( mod(n,2)~=0 );
  h(oddNs) = -1 ./ ( n(oddNs).*n(oddNs) * pi*pi * dSize*dSize );
  zeroIndx = find( n==0 );
  h(zeroIndx) = 1/(4*dSize*dSize);
  nPadded = 2*nDetectors;
  hZP = zeros(1,nPadded);
  hZP(1:nDetectors) = h;

  if strcmp( window, 'Hanning' )
    df = 1 / ( dSize * nPadded );
    halfP = floor(nPadded / 2) + 1;
    fIndxs = ( (1:nPadded) - halfP );
    f = df .* fIndxs;
    absF = abs( f );

    fN = 1 ./ ( 2 * dSize );
    hannWin = 0.5 .* ( 1 + cos( pi .* absF ./ fN ) );
    hannWin = hannWin .* ( absF < fN );
    fftHZP = fftshift( fft(hZP) );
    filtFftHZP = fftHZP .* hannWin;
    hZP =  ifft( ifftshift( filtFftHZP) );
  end

  hZP = ones( nThetas, 1 ) * hZP;
  sinoZP = zeros( nThetas, nPadded );
  sinoZP( :, 1:nDetectors ) = sinogram;

  fftH = fft( hZP, [], 2 );
  fftSino = fft( sinoZP, [], 2 );

  fftFiltSino = fftH .* fftSino;
  filtSino = dSize * ifft( fftFiltSino, [], 2 );
  filtSino = circshift( filtSino, [0 -zeroIndx+1] );
  filtSino = filtSino( :, 1:nDetectors );

  %% Perform backprojection
  out = ctBackProject( filtSino, thetas, dSize, cx, cy, Nx, Ny, dx, dy );
  out = out * pi / nThetas;
end
