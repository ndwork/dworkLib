
function out = projectVolume( vol, H, varargin )
  % out = projectImageWithHomography( img, H, [ range ] )
  %
  % Inputs:
  % vol - a 3D array of data to project
  % H - the 4x4 homography
  % range - An optional 4D array specifying the range to project to
  %   [xmin xmax ymin ymax zmin zmax]
  %
  % Written by Nicholas Dwork 2016

  sVol = size( vol );

  defaultRange = [ 1 sVol(2) 1 sVol(1) 1 sVol(3) ];
  p = inputParser;
  p.addOptional( 'range', defaultRange );
  p.parse( varargin{:} );
  range = p.Results.range;

  N = numel( vol );
  xs = range(1):range(2);
  ys = range(3):range(4);
  zs = range(5):range(6);
  [xs, ys, zs] = meshgrid( xs, ys, zs );
  coords = [ xs(:)'; ys(:)'; zs(:)'; ones(1,N) ];

  projCoords_h = inv(H) * coords;
  projCoords = hom2Euc( projCoords_h );
  interped = interp3( xs, ys, zs, vol, ...
    projCoords(1,:), projCoords(2,:), projCoords(3,:), ...
    'linear', 0);
  out = reshape( interped, sVol );

end
