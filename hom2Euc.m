
function euc = hom2Euc( hom )
  % euc = hom2Euc( hom )
  % Convert homogeneous coordinates to Euclidean coordinates
  %
  % Inputs:
  % hom - 2D array of size (N+1)xM where N is the number of dimensions
  %   representing the homogeneous coordinates of a set of points,
  %   and M is the number of points.
  %
  % Outputs:
  % euc - 2D array of size NxM
  %
  % Written by Nicholas Dwork 2016

  sHom = size(hom);

  euc = zeros( sHom(1)-1, sHom(2) );

  for i=1:sHom(1)-1
    euc(i,:) = hom(i,:) ./ hom( sHom(1), : );
  end

end
