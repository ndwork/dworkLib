
function out = volMatrixProd( vol, A )
  % out = volMatrixProd( vol, A )
  %
  % Inputs:
  % vol - 3D array of size NxPxK
  % M - 2D array of size MxN
  %
  % Left matrix multiply A by each slice of vol
  %
  % Written by Nicholas Dwork - Copyright 2018
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp( 'Usage:  out = volMatrixProd( vol, A )' );
    return;
  end

  sVol = ones(1,3);
  sVol(1:ndims(vol)) = size( vol );
  out = zeros( [ sVol(1) size(A,2) sVol(3) ] );

  if sVol(3) < 2 * sVol(1)*sVol(2)
    for m=1:sVol(3)
      out(:,:,m) = vol(:,:,m) * A;
    end

  else

    volT = permute( vol, [2 1 3] );
    rVolT = reshape( volT, [ sVol(2), sVol(1)*sVol(3) ] );

    for u=1:size(A,2)
      tmp = sum( bsxfun( @times, A(:,u), rVolT ), 1 );
      out(:,u,:) = reshape( tmp, [sVol(1) 1 sVol(3)] );
    end
  end
end

