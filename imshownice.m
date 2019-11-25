
function [scaling,imH] = imshownice( img, varargin )
  % scaling = imshownice( img [, scale, 'method', method, 'sdevScale', sdevScale, ...
  %   'border', border ] )
  % show the image on the following scale:
  %   meanImg - sdevScale*sdevImg, meanImg + sdevScale*sdevImg
  %
  % Inputs:
  % img - 2D array representing the image
  %
  % Optional Inputs:
  % scale - factor to scale the size of the image for display
  % method - when scaling method, interpolation method to use
  %   default is 'nearest'
  %   any method accepted by imresize is accepted for this parameter
  % sdevScale - default is 2.5
  % border - border to put around the image in the figure window
  %   either 'noBorder' or a value in pixels  (default is 10)
  %
  % Outputs:
  % scaling - the range of intensities used in the display.
  % imH - image handle
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  defaultScale = 1.0;
  defaultMethod = 'nearest';
  defaultSDevScale = 2.5;
  defaultBorder = 10;
  p = inputParser;
  p.addOptional( 'scale', defaultScale, @isnumeric );
  p.addParameter( 'method', defaultMethod, @(x) true );
  p.addParameter( 'sdevScale', defaultSDevScale, @isnumeric );
  p.addParameter( 'border', defaultBorder );
  p.parse( varargin{:} );
  scale = p.Results.scale;
  method = p.Results.method;
  sdevScale = p.Results.sdevScale;
  border = p.Results.border;

  if numel( sdevScale ) == 0, sdevScale = defaultSDevScale; end
  
  medianImg = median( real(img(:)) );
  sdevImg = std( real(img(:)) );

  if ismatrix( img )
    % Grayscale image
    tmp = imresize( img, scale, method );
    scaling = [ medianImg - sdevScale*sdevImg, medianImg + sdevScale*sdevImg ];
    imH = imshow( tmp, scaling );
  else
    % Color image
    sImg = size(img);
    tmp = zeros( [ sImg(1:end-1) * scale, sImg(end) ] );
    for i=1:sImg(end)
      tmp(:,:,i) = imresize( img(:,:,i), scale, method );
    end
    inMin = medianImg - sdevScale*sdevImg;
    inMax = medianImg + sdevScale*sdevImg;
    scaling = [inMin inMax];
    scaled = scaleImg( tmp, scaling, [0 1] );
    imH = imshow( scaled );
  end

  if ischar( class(border) ) && strcmp( border, 'none' )
    displayBorder = 0;
  elseif border < 0
    error( 'border must be great than or equal to 0.' );
  else
    displayBorder = border;
  end
  ca = gca;
  beforeAxesUnits = ca.Units;
  set( ca, 'units', 'pixels' );
  x = get( ca, 'position' );
  cf = gcf;
  beforeFigUnits = cf.Units;
  set( cf, 'units', 'pixels' );
  y = get( cf, 'position' );
  if( displayBorder > 0 )
    set( cf, 'position', [y(1) y(2) x(3)+2*displayBorder x(4)+2*displayBorder+10] );
      % set the position of the figure to the length and width of the axes
      % add 10 pixels to make space for a title
    set( ca, 'position', [displayBorder displayBorder x(3) x(4)] );
  else
    set( cf, 'position', [y(1) y(2) x(3) x(4)] );
    set( ca, 'position', [0 0 x(3) x(4)] );
  end
  % Now restore units to previously used values
  set( ca, 'units', beforeAxesUnits );
  set( cf, 'units', beforeFigUnits );

  drawnow; pause(0.02);
  
  addToolbarExplorationButtons( cf );  % Restore the missing toolbar
end
