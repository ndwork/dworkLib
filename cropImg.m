
function out = cropImg( img, N )
  % out = cropImg( img, N )
  % Crops out the center region of the image.
  % (0,0) is defined according to fftshift
  %
  % Inputs:
  % img - 2D array
  % N - specified the size of the cropped image
  %   If N is a scalar, then a square image is extracted
  %   If N is a 2 element array, then the final size is [N(1) N(2)]

  if numel(N)==2
    Ny=N(1); Nx=N(2);
  else
    Ny=N; Nx=N;
  end

  [My, Mx] = size(img);

  % Find y portion
  halfMy = My/2;
  if mod(My,2)==0
    cy = halfMy + 1;
    if mod(Ny,2)==0
      halfNy = Ny/2;
      minY = cy - halfNy;
      maxY = cy + halfNy - 1;
    else
      halfNy = floor(Ny/2);
      minY = cy - halfNy;
      maxY = cy + halfNy;
    end
  else
    cy = ceil(halfMy);
    if mod(Ny,2)==0
      halfNy = Ny/2;
      minY = cy - halfNy;
      maxY = cy + halfNy - 1;
    else
      halfNy = floor(Ny/2);
      minY = cy - halfNy;
      maxY = cy + halfNy;
    end
  end

  % Find x portion
  halfMx = Mx/2;
  if mod(Mx,2)==0
    cx = halfMx + 1;
    if mod(Nx,2)==0
      halfNx = Nx/2;
      minX = cx - halfNx;
      maxX = cx + halfNx - 1;
    else
      halfNx = floor(Nx/2);
      minX = cx - halfNx;
      maxX = cx + halfNx;
    end
  else
    cx = ceil(halfMx);
    if mod(Nx,2)==0
      halfNx = Nx/2;
      minX = cx - halfNx;
      maxX = cx + halfNx - 1;
    else
      halfNx = floor(Nx/2);
      minX = cx - halfNx;
      maxX = cx + halfNx;
    end
  end
  
  out = img( minY:maxY, minX:maxX );
end
