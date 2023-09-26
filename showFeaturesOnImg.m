
function showFeaturesOnImg( features, varargin )
  % showFeaturesOnImg( features [, img, 'range', range, 'scale', scale, ...
  %   'color', color, 'figH', figH, 'connect', true/false, 'offset', offset ] )
  %
  % Inputs:
  % features - 2D array of size Nx2
  %   N is the number of points
  %   The first/second column is the x/y location
  %
  % Optional Inputs:
  % img - if included, also displays the image
  % range - 2 element array specifying image display range ( default is [] )
  % scale - magnify image by this scale
  % connect - if set to true, draws lines between features
  % offset - 2 element array specifying the amount to add to each feature to
  %   convert to image coordinates for plot
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addOptional( 'img', [] );
  p.addParameter( 'connect', false, @(x) islogical(x) );
  p.addParameter( 'figH', [] );
  p.addParameter( 'LineWidth', 2 );
  p.addParameter( 'MarkerSize', 5 );
  p.addParameter( 'offset', [], @(x) isnumeric(x) && numel(x) == 2 );
  p.addParameter( 'range', [] );
  p.addParameter( 'scale', 1 );
  p.addParameter( 'color', 'y' );
  p.parse( varargin{:} );
  img = p.Results.img;
  connect = p.Results.connect;
  figH = p.Results.figH;
  LineWidth = p.Results.LineWidth;
  MarkerSize = p.Results.MarkerSize;
  offset = p.Results.offset;
  scale = p.Results.scale;
  range = p.Results.range;
  color = p.Results.color;

  if numel( offset ) == 0
    offset = -0.5 * ones( size( features ) );
  else
    offset = repmat( offset, [size(features,1) 1] );
  end
  
  if numel( img ) > 0
    if numel( figH ) > 0
      figure( figH );
    else
      figH = figure;
    end
    imshowscale( img, scale, 'range', range );
  end
  if numel( figH ) > 0, figure( figH ); end
  hold on;
  
  rFeatures = round( scale * ( features + offset ) );
  plot( rFeatures(:,1), rFeatures(:,2), [color,'x'], 'LineWidth', LineWidth, 'MarkerSize', MarkerSize );
  drawnow;
  
  if connect == true
    hold on;
    plotnice( rFeatures(:,1), rFeatures(:,2) );
  end
end
