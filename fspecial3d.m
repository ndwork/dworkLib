
function h = fspecial3d( type, varargin )
  % out = fspecial3d( type, parameters )
  %
  % Creates a 3D filter h of the specified type.  The different types are
  % listed below.
  %
  % h = fspecial3d( 'average', hSize )
  %   hSize can be a 3 element array or a scalar
  %
  % Written by Nicholas Dwork

  switch type
    case 'average'
      h = fspecial3d_average( varargin{:} );
    case 'gaussian'
      h = fspecial3d_gaussian( varargin{:} );
  end
end

function h = fspecial3d_average( hSize )
  if numel(hSize)==1
    ny=hSize; nx=hSize; nz=hSize;
  else
    ny=hSize(1); nx=hSize(2); nz=hSize(3);
  end
  N = ny * nx * nz;
  h = 1/N * ones( ny, nx, nz );
end

function h = fspecial3d_gaussian( hSize, sigma )
  if numel(hSize)==1
    ny=hSize; nx=hSize; nz=hSize;
  else
    ny=hSize(1); nx=hSize(2); nz=hSize(3);
  end

  halfY = floor(ny/2);
  if mod(ny,2)==0
    y = -halfY:halfY-1;
  else
    y = -halfY:halfY;
  end

  halfX = floor(nx/2);
  if mod(nx,2)==0
    x = -halfX:halfX-1;
  else
    x = -halfX:halfX;
  end

  halfZ = floor(nz/2);
  if mod(nz,2)==0
    z = -halfZ:halfZ-1;
  else
    z = -halfZ:halfZ;
  end

  [x,y,z] = meshgrid(x,y,z);
  h = exp( -0.5/(sigma*sigma) * ( x.*x + y.*y + z.*z ) );
  h = h / sum(h(:));
end

