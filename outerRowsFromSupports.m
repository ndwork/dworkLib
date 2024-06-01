
function [ outerRows, innerMask ] = outerRowsFromSupports( supports )
  % [ outerRows, innerMask ] = outerRowsFromSupports( supports )
  %
  % This is a supporting function for reconNonRectSupports.
  %
  % Written by Nicholas Dwork, Copyright 2024
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  nRows = size( supports, 1 );
  nCols = size( supports, 2 );
  shiftedSupports = circshift( supports, floor(nCols/2), 2 );
  summedSupports = supports + shiftedSupports;

  if nargout > 1
    innerMask = ( summedSupports == 1 ) .* supports;
  end

  nCoils = size( supports, 3 );
  outerRows = zeros( nRows, nCoils );

  for coil = 1 : nCoils
    summedSupport = summedSupports(:,:,coil);
    if max( summedSupport(:) ) == 1
      outerRows(:,coil) = 1;
      continue;
    end

    maxRowsSummedSupport = max( summedSupports(:,:,coil), [], 2 );
    firstOuterRowIndx = find( maxRowsSummedSupport > 1, 1, 'first' ) - 1;
    lastOuterRowIndx = find( maxRowsSummedSupport > 1, 1, 'last' ) + 1;
    
    outerRows( 1 : firstOuterRowIndx, coil ) = 1;
    outerRows( lastOuterRowIndx : end, coil ) = 1;
  end
end
