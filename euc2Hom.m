
function hom = euc2Hom( euc )
  % hom = hom2Euc( euc )
  % Convert Euclidean coordinates to homogeneous coordinates
  %
  % Inputs:
  % euc - 2D array of size NxM where M is the number of points and N is
  %   the number dminensions.
  %
  % Outputs:
  % hom - 2D array of size (N+1)xM.
  %
  % Written by Nicholas Dwork 2016

  sEuc = size(euc);
  
  hom = ones( sEuc(1)+1, sEuc(2) );
  hom(1:sEuc(1),:) = euc;
end