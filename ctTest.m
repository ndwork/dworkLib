

function ctTest
  % Demo routine for ct Codes
  %
  % Written by Nicholas Dwork and Uzair Sikora, Copyright 2013
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  clear; close all; rng(1);

  im = phantom();
  delta = 0.001;
  nDetectors = 500;
  dSize = 0.001;
  dTheta = 1 * pi/180;
  thetas = 0:dTheta:pi-dTheta;

  % Reconstruction parameters
  cx = 0;   Nx=256;   dx=delta;
  cy = 0;   Ny=256;   dy=delta;


  figure;  imshowscale( im );  titlenice( 'Original Image' );

  sinogram = ctRadon( im, delta, nDetectors, dSize, thetas );
  figure;  imshowscale( sinogram );  titlenice( 'sinogram' );

  bp = ctBackProject( sinogram, thetas, dSize, cx, cy, Nx, Ny, dx, dy );
  figure;  imshowscale( bp );  titlenice( 'Back Projection' );

  recon = ctIRadon( sinogram, thetas, dSize, cx, cy, Nx, Ny, dx, dy, 'Hanning' );
  figure;  imshowscale( recon );  titlenice( 'Reconstruction' );

  recon2DFT = ctIRadon2DFT( sinogram, thetas, dSize, cx, cy, Nx, Ny, dx, dy, 'none' );
  figure;  imshowscale( recon2DFT );  titlenice( '2DFT Reconstruction' );

  sinogram2 = ctRadon( recon, delta, nDetectors, dSize, thetas );
  figure;  imshowscale( sinogram2 );  titlenice( 'Singoram 2' );

  u = rand( Ny, Ny );
  Ru = ctRadon( u, delta, nDetectors, dSize, thetas );
  v = rand( size(Ru) );
  Ru_v = (Ru(:))' * v(:);
  Rtv = ctBackProject( v, thetas, dSize, cx, cy, Nx, Ny, dx, dy );
  u_Rtv = u(:)' * Rtv(:);
  diff = Ru_v - u_Rtv;
  disp(['Difference is: ', num2str(diff)]);

end
