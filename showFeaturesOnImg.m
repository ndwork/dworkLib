
function showFeaturesOnImg( features, varargin )
  % showFeaturesOnImg( features [, img, 'range', range, 'scale', scale, 'color', color] )
  %
  % Inputs:
  % features - 2D array of size Nx2
  %   N is the number of points
  %   The first/second column is the x/y location
  %
  % Optional Inputs:
  %   img - if included, also displays the image
  %   range - 2 element array specifying image display range ( default is [] )
  %   scale - magnify image by this scale
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addOptional( 'img', [] );
  p.addParameter( 'range', [] );
  p.addParameter( 'scale', 1 );
  p.addParameter( 'color', 'y' );
  p.parse( varargin{:} );
  img = p.Results.img;
  scale = p.Results.scale;
  range = p.Results.range;
  color = p.Results.color;

  if numel( img ) > 0
    figure; imshowscale( img, scale, 'range', range );
  end
  hold on
  rFeatures = round( scale * features );
  plot( rFeatures(:,1), rFeatures(:,2), [color,'*']);
  drawnow;
end
