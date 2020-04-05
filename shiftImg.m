
function out = shiftImg( img, shifts )
  % out = shiftImg( img, shifts )
  %
  % Shift the image; shifts are defined according to circshift
  % Regions with unknown data are zero filled.
  %
  % Inputs:
  % img - a 2D array
  % shifts - a 2 element array specifying vertical and horiztonal shift
  %
  % Outputs:
  % out - the shifted image
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  out = circshift( img, shifts );

  if shifts(1) > 0
    out(1:shifts(1),:) = 0;
  elseif shifts(1) < 0
    out(end+shifts(1)+1:end,:) = 0;
  end

  if shifts(2) > 0
    out(:,1:shifts(2)) = 0;
  elseif shifts(2) < 0
    out(:,end+shifts(2)+1:end) = 0;
  end

end

