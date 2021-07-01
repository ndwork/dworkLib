
function contour = drawContour( img, varargin )
  % contour = drawContour( img )
  %
  % displays figure to the screen and lets user create a contour
  %   Type 'h' for help
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

  p = inputParser;
  p.addOptional( 'scale', 1, @isnumeric );
  p.parse( varargin{:} );
  scale = p.Results.scale;

  figure;  imshowscale( img, scale );
  hold on;  grid on;

  contour = [];
  while true
    title( 'type h for help' );
    try
      [y,x,button] = ginput(1);
      zoom out;
    catch
      button = [];
    end

    if numel( button ) == 0
      % Figure was closed before key was pressed
      break;

    elseif button == 'h'
      disp( 'Left click with mouse to select a contour point.');
      disp( 'h for this help' );
      disp( 'q to quit the application' );
      disp( 'z to zoom in and out of the image' );

    elseif button == 'q'
      close( gcf );
      break;

    elseif button == 'z'
      zoom on;
      title( 'Hit any key to end zoom' );
      pause();  zoom off;

    else
      contour = [ contour; [ y x ]; ];   %#ok<AGROW>
      plot( contour(:,1), contour(:,2), 'yx' );
      plot( contour(:,1), contour(:,2), 'y' );
    end

  end

  if scale ~= 1, contour = contour / scale; end
end
