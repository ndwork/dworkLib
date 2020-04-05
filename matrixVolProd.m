
function out = matrixVolProd( A, vol )
  % out = matrixVolProd( A, vol )
  %
  % Left matrix multiply each slice of vol by the matrix A.
  % If the 3rd dimension of vol is much larger than the dimensions of A, this
  % is faster than just writing a loop.
  %
  % Inputs:
  % A - 2D array of size MxN
  % vol - 3D array of size NxPxK
  %
  % Written by Nicholas Dwork - Copyright 2018
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 2
    disp( 'Usage:  out = matrixVolProd( A, vol )' );
    return;
  end
  
  sVol = ones(1,3);
  sVol(1:ndims(vol)) = size( vol );
  out = zeros( [ size(A,1) sVol(2) sVol(3) ] );

  if sVol(3) < 2 * sVol(1)*sVol(2)

    for m=1:sVol(3)
      out(:,:,m) = A * vol(:,:,m);
    end

  else

    AT = transpose( A );
    rVol = reshape( vol, [ sVol(1), sVol(2)*sVol(3) ] );

    for u=1:size(A,1)
      tmp = sum( bsxfun( @times, AT(:,u), rVol ), 1 );
      out(u,:,:) = reshape( tmp, [1 sVol(2) sVol(3)] );
    end
  end
end

