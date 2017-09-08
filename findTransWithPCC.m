
function [vShift,hShift] = findTransWithPCC( img1, img2 )
  % [vShift,hShift] = findTransWithPCC( img1, img2 )
  %
  % Outputs:
  % vShift - the vertical shift to apply to img1 so that it is aligned
  %          with image 2
  % hShift - the horizontal shift to apply to img1 so that it is aligned
  %          with image 2
  %
  % Written by Nicholas Dwork - Copyright 2017
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  pcc = phaseCrossCorrelate( img1, img2 );
  [~,maxIndx] = max( pcc(:) );
  [vShift,hShift] = ind2sub( size(pcc), maxIndx );
  vShift = vShift - 1;
  hShift = hShift - 1;
end