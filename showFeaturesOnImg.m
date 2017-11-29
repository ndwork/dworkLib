
function showFeaturesOnImg( features, img, varargin )
  % showFeaturesOnImg( features, img [, range, 'scale', scale, 'color', color] )
  %
  % Inputs:
  % features - 2D array of size Nx2
  %   N is the number of points
  %   The first/second column is the x/y location
  %
  % Optional Inputs:
  %   range - 2 element array specifying image display range
  %     ( default is [0 1] )
  %   scale - magnify image by this scale
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addOptional( 'range', [0 1] );
  p.addParameter( 'scale', 1 );
  p.addParameter( 'color', 'y' );
  p.parse( varargin{:} );
  scale = p.Results.scale;
  range = p.Results.range;
  color = p.Results.color;

  if nSigDims( features ) > 1
    theseFeatures = features;
  else
    [theseFeaturesY,theseFeaturesX] = ind2sub( size(img), features );
    theseFeatures = [ theseFeaturesY', theseFeaturesX' ];
  end

  figure;
  imshowscale( img, scale, 'range', range );
  hold on
  rFeatures = round( scale * theseFeatures );
  plot( rFeatures(:,1), rFeatures(:,2), [color,'*']);
  drawnow;
end
