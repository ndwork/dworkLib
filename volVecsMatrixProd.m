
function out = volVecsMatrixProd( A, vs )
  % out = volVecsMatrixProd( A, vs )
  %
  % Matrix multiplies each slice of A by each colum of vs
  %
  % Inputs:
  % A - A 3D array
  % vs - A 2D array representing a set of vectors (each col is a vector)
  %
  % Outputs:
  % out - out(:,i) = A(:,:,i) * vs(:,i) for all i
  %
  % Written by Nicholas Dwork - Copyright 2018
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  out = zeros( size(A,1), size(vs,2) );
  for i=1:size(A,2)
    out = out + squeeze(A(:,i,:)) .* repmat( vs(i,:), [size(A,1) 1] );
  end

end

