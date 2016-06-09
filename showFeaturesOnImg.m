
function showFeaturesOnImg( features, img, varargin )
  % showFeaturesOnImg( features, img [, scale] )
  %
  % Inputs:
  % features - 2D array of size Nx2
  %   N is the number of points
  %   The first/second column is the x/y location
  %
  % Optional Inputs:
  %   scale - 2 element array specifying image scale
  %
  % Written by Nicholas Dwork - Copyright 2016

  p = inputParser;
  p.addOptional( 'scale', [] );
  p.parse( varargin{:} );
  scale = p.Results.scale;
  
  figure;
  imshow(img, scale);
  hold on
  plot( features(:,1), features(:,2), 'r*');
  drawnow;
end
