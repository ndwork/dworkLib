

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
  %   horiontal directions that must be applied to image 1 to align it with
  %   image 2
  %
  % Written by Nicholas Dwork - Copyright 2017
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  pc = phaseCrossCorrelate( img2, img1 );
  [~,maxIndx] = max(pc(:));
  sImg = size( img1 );
  [shiftedY, shiftedX] = ind2sub(sImg, maxIndx );
  shiftedY = shiftedY - 1;  % Account for 1 based indexing of Matlab
  shiftedX = shiftedX - 1;
  shifts = [ shiftedY shiftedX ];
end

