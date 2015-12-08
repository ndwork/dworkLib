
function [ky,kx] = size2fftCoordinates( N )
  % [ky,kx] = size2fftCoordinates( N )
  %
  % Inputs:
  %   N is a two element array [Ny Nx] specifying the number of row and
  %     columns of the Fourier Transformed image
  %
  % Outputs:
  %   Let fftImg = fftshift( fft2( ifftshift( img ) ) );
  %     ky is a vector of size( img, 1 ) with the k-space location of
  %       each column of fftImg
  %     kx is a vector of size( img, 1 ) with the k-space location of
  %       each row of fftImg

  
  
  Ny = N(1);
  dky = 1 / Ny;
  if mod( Ny, 2 ) == 0
    ky = linspace(-0.5,0.5-dky,Ny);
  else
    ky = linspace(-0.5+dky/2,0.5-dky/2,Ny);
  end

  Nx = N(2);
  dkx = 1 / Nx;
  if mod( Nx, 2 ) == 0
    kx = linspace(-0.5,0.5-dkx,Nx);
  else
    kx = linspace(-0.5+dkx/2,0.5-dkx/2,Nx);
  end
end