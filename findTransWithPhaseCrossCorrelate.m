

function shifts = findTransWithPhaseCrossCorrelate( img1, img2 )
  % [shiftedY, shiftedX] = findTransWithPhaseCrossCorrelate( img1, img2 )
  % Finds the translations between images using phase cross correlation
  %
  % Inputs:
  % img1 - a 2D array representing an image
  % img2 - a 2D array representing an image
  %
  % Output:
  % shifts = A two element array specifying the shifts in the vertical and
  %   horiontal directions
  %
  % Written by Nicholas Dwork - Copyright 2017
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  pc = phaseCrossCorrelate( img1, img2 );
  [~,maxIndx] = max(pc(:));
  sImg = size( img1 );
  [shiftedY, shiftedX] = ind2sub(sImg, maxIndx );
  shiftedY = mod( -(shiftedY-1), sImg(1) );
  shiftedX = mod( -(shiftedX-1), sImg(2) );
  shifts = [ shiftedY shiftedX ];
end

