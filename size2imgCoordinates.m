
function [y,x] = size2imgCoordinates( N )
  % [ky,kx] = size2imgCoordinates( N )
  % Determine the image coordinates where (0,0) is at the "center"
  %   (center is defined in the same way as fftshift)
  %
  % Inputs:
  %   N is a two element array [Ny Nx] specifying the number of row and
  %     columns of the image
  %
  % Outputs:
  %   y is a vector representing the y locations of each row of the image
  %   x is a vector representing the x locaitons of each row of the image


  Ny = N(1);
  y = (0:Ny-1) - floor(0.5*Ny);

  Nx = N(2);
  x = (0:Nx-1) - floor(0.5*Nx);
end
