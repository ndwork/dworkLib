
function imH = imshowscale( img, varargin )
  % imshowscale( img, [ scale, , 'border', border, 'method', method, 'range', range, ...
  %   'thresh', thresh ] )
  % displays figure to the screen where size of image is scaled by scale
  %
  % Inputs:
  % img - 2D array representing the grayscale image or 3D array
  %       representing the color image (third dimension has size 3)
  %
  % Outputs:
  % imH - a handle to the image object
  %
  % Optional Inputs:
  % scale - factor to scale the size of the image for display (default is 1)
  % method - when scaling method, interpolation method to use
  %   default is 'nearest'
  %   any method accepted by imresize is accepted for this parameter
  % range - two element array specifying the display range of intensities
  %   or an image
  %   If range is [], sets equal to [min(img(:)) max(img(:))]
  %     This is the default.
  %   If range is 'nice', uses imshownice to display image
  %   If range is an image (larger than two elements), then the range of
  %     the display equals the dynamic range of the input array
  % border - border to put around the image in the figure window
  %   either 'noBorder' or a value in pixels (default is 10) or an array
  %   with four elements [ left bottom right top ]
  %
  % Written by Nicholas Dwork - Copyright 2017
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  defaultScale = 1;
  defaultMethod = 'nearest';
  defaultBorder = 10;
  defaultFontSize = 20;
  p = inputParser;
  p.addOptional( 'scale', defaultScale, @isnumeric );
  p.addParameter( 'border', defaultBorder );
  p.addParameter( 'fontSize', defaultFontSize, @isnumeric );
  p.addParameter( 'method', defaultMethod );
  p.addParameter( 'range', [] );
  p.addParameter( 'sdevScale', [], @isnumeric );
  p.addParameter( 'thresh', [], @ispositive );
  p.addParameter( 'xAxis', [], @isnumeric );
  p.addParameter( 'yAxis', [], @isnumeric );
  p.parse( varargin{:} );
  scale = p.Results.scale;
  border = p.Results.border;
  fontSize = p.Results.fontSize;
  method = p.Results.method;
  range = p.Results.range;
  sdevScale = p.Results.sdevScale;
  thresh = p.Results.thresh;
  xAxis = p.Results.xAxis;
  yAxis = p.Results.yAxis;

  img = double(img);

  if figureExists()
    cf = gcf();
    
    thisFrame = getframe( cf );
    thisFrameData = thisFrame.cdata;
    
    if min( thisFrameData(:) ) ~= 240 || max( thisFrameData(:) ) ~= 240
      outerPos = cf.OuterPosition;
    end
  end

  if strcmp( 'nice', range )
    [~,imH] = imshownice( img, scale, 'method', method, 'sdevScale', sdevScale, 'thresh', thresh );

  else

    if ismatrix( img )
      if numel( range ) == 0
        thisRange = range;
      else
        thisRange = [ min(range(:)) max(range(:)) ];
      end

      imH = imshow( imresize( img, scale, method ), thisRange );

    elseif ndims(img) == 3

      if numel( range ) > 0
        minRange = min( range );
        maxRange = max( range );
        dRange = maxRange - minRange;
        img = ( img - minRange ) / dRange;
      end

      colorResized = imColorResize( img, scale, method );
      scaled = scaleImg( colorResized, [0 1], [ min(img(:)) max(img(:)) ] );
      imH = imshow( scaled );

    else
      error('wrong number of dimensions of img');

    end
  end


  if ischar( class(border) ) && strcmp(border,'none')
    displayBorder = zeros(1,4);
  elseif border < 0
    error('border must be great than or equal to 0.');
  elseif numel( border ) == 1
    displayBorder = border * ones(1,4);
  elseif numel( border ) == 4
    displayBorder = border;
  end
  ca = gca;
  beforeAxesUnits = ca.Units;
  set( ca, 'units', 'pixels' );
  x = get( ca, 'position' );
  cf = gcf;
  beforeFigUnits = cf.Units;
  set( cf, 'units', 'pixels' );
  y = get(cf,'position');
  if( numel(displayBorder) > 0 || displayBorder > 0 )
    set( cf, 'position', [y(1) y(2) ...
      x(3)+displayBorder(1)+displayBorder(3) ...
      x(4)+displayBorder(2)+displayBorder(4)+10] );
      % set the position of the figure to the length and width of the axes
      % add 10 pixels to make space for a title
    set( ca, 'position', [displayBorder(1) displayBorder(2) x(3) x(4)]);
  else
    set( cf,'position', [y(1) y(2) x(3)-1 x(4)-1] );
    set( ca,'position', [0 0 x(3) x(4)]);
  end

  % If desired, draw axes
  if numel( xAxis ) + numel( yAxis ) > 0
    xTextBuffer = ( numel( xAxis ) > 0 ) * 2*fontSize + 10;
    yTextBuffer = ( numel( yAxis ) > 0 ) * 2*fontSize;
    axis on;
    set( cf, 'units', 'pixels' );
    set( ca, 'units', 'pixels', 'fontsize', fontSize );
    cfPos = get( cf, 'position' );    % Position is [ left bottom width height ]
    caPos = get( ca, 'position' );
    newCfPos = [  cfPos(1), cfPos(2), cfPos(3)+xTextBuffer, cfPos(4)+yTextBuffer ];
    newCaPos = [  caPos(1)+xTextBuffer, caPos(2)+yTextBuffer, caPos(3), caPos(4) ];
    set( cf, 'position', newCfPos );
    set( ca, 'position', newCaPos, 'box', 'off' );

    set( cf, 'units', beforeFigUnits );
    set( ca, 'units', beforeAxesUnits );
  end

  if exist( 'outerPos', 'var' )
    thisOuterPos = cf.OuterPosition;
    thisOuterPos(1:2) = outerPos(1:2);
    cf.OuterPosition = thisOuterPos;
  end
  
  drawnow; pause(0.02);

  % Now restore units to previously used values
  set( ca, 'units', beforeAxesUnits );
  set( cf, 'units', beforeFigUnits );

  addToolbarExplorationButtons( cf );  % Restore the missing toolbar.
end
