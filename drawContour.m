function contour = drawContour( img )
  % contour = drawContour( img )
  % displays figure to the screen and lets user create a contour
  %   Left click with the mouse to select a new point on the contour
  %   Hit 'q' to complete contour selection
  %
  % Inputs:
  % img - 2D array representing the grayscale image or 3D array
  %       representing the color image (third dimension has size 3)
  %
  % Outputs:
  % contour - 2D array of size N x 2 where N is the number of contour
  %   points
  %
  % Written by Jolie Wang - Copyright 2021
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp( 'Usage:  contour = drawContour( img )' );
    contour = [];
    return;
  end

  figure;  imshow( img );  hold on;  grid on;
  disp('Press q to quit or mouse click to select contour points: ');

  contour = [];
  while true
    [y,x,button] = ginput(1);

    if (button == double('q')), break; end
    
    contour = [ contour; [ y x ]; ];   %#ok<AGROW>

    plot( contour(:,1), contour(:,2), 'yx' );
    plot( contour(:,1), contour(:,2), 'y' );
  end

  close( gcf );
end
