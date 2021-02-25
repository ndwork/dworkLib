
% Demo routine for ct Codes
% Written by Nicholas Dwork and Uzair Sikora


function ctTest
  clear; close all;

  im = phantom();
  delta = 0.001;

  figure( 'name', 'Original Image' );
  imshow( im, [] );

  
  nDetectors = 500;
  dSize = 0.001;
  dTheta = 1 * pi/180;
  thetas = 0:dTheta:pi-dTheta;

  % Reconstruction parameters
  cx = 0;   Nx=256;   dx=delta;
  cy = 0;   Ny=256;   dy=delta;


  type = 'fast';
  sinogram = ctRadon( im, delta, nDetectors, dSize, thetas, type );

  figure( 'name', 'sinogram' );
  imshow( sinogram, [] );


  bp = ctBackProject( sinogram, thetas, dSize, cx, cy, Nx, Ny, ...
    dx, dy, type );
  figure( 'name', 'Back Projection' );
  imshow( bp, [] );

  recon = ctIRadon( sinogram, thetas, dSize, cx, cy, Nx, Ny, dx, dy, 'Hanning' );
  figure( 'name', 'Reconstruction' );
  imshow( recon, [] );

  recon2DFT = ctIRadon2DFT( sinogram, thetas, dSize, cx, cy, Nx, Ny, dx, dy, 'none' );
  figure( 'name', '2DFT Reconstruction' );
  imshow( recon2DFT, [] );

  sinogram2 = ctRadon( recon, delta, nDetectors, dSize, thetas );
  figure( 'name', 'Sinogram 2' );
  imshow( sinogram2, [] );
  
  
  type = 'iso';
  u = rand( Ny, Ny );
  Ru = ctRadon( u, delta, nDetectors, dSize, thetas, type );
  v = rand( size(Ru) );
  Ru_v = (Ru(:))' * v(:);
  Rtv = ctBackProject( v, thetas, dSize, cx, cy, Nx, Ny, ...
    dx, dy, type );
  u_Rtv = u(:)' * Rtv(:);
  diff = Ru_v - u_Rtv;
  disp(['Difference is: ', num2str(diff)]);
  
  
end
