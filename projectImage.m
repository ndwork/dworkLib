
function out = projectImage( img, H, varargin )
  % out = projectImageWithHomography( img, H, [ range ] )
  %
  % Inputs:
  % img - a 2D array of data to project
  % H - the 3x3 homography
  % range - An optional 4D array specifying the range to project to
  %   [xmin xmax ymin ymax]
  %
  % Written by Nicholas Dwork 2016

  sImg = size( img );

  defaultRange = [ 1 sImg(2) 1 sImg(1) ];
  p = inputParser;
  p.addOptional( 'range', defaultRange );
  p.parse( varargin{:} );
  range = p.Results.range;

  N = numel( img );
  xs = range(1):range(2);
  ys = range(3):range(4);
  [xs, ys] = meshgrid( xs, ys );
  coords = [ xs(:)'; ys(:)'; ones(1,N) ];

  projCoords_h = inv(H) * coords;
  projCoords = hom2Euc( projCoords_h );
  interped = interp2( xs, ys, img, projCoords(1,:), projCoords(2,:), ...
    'linear', 0);
  out = reshape( interped, sImg );

end
