
function out = rotateAndTranslateImg( R, t, img, varargin )
  % out = rotateAndTranslateImg( R, t, img [, extrapVal] )
  %
  % Inputs:
  % R - a 2x2 rotation matrix
  % t - a two element array specifying translation in pixels
  % img - a 2D array
  %
  % Optional Inputs:
  % extrapVal - value to put where extrapolation takes place
  %
  % Outputs:
  % out - a 2D array - the transformed image
  %
  % Rotates and translates image according to the formula
  %   newPt = R * pt + t;
  %
  % Written by Nicholas Dwork - Copyright 2016

  defaultExtrapVal = 0;
  p = inputParser;
  p.addOptional( 'extrapVal', defaultExtrapVal );
  p.parse( varargin{:} );
  extrapVal = p.Results.extrapVal;

  sImg = size( img );
  ys = 1:sImg(1);
  xs = 1:sImg(2);
  [xs, ys] = meshgrid( xs, ys );
  coords = [ xs(:)'; ys(:)' ];

  tCoords = zeros( size(coords) );
  tCoords(1,:) = coords(1,:) - t(1);
  tCoords(2,:) = coords(2,:) - t(2);
  RCoords = inv(R) * tCoords;
  interped = interp2( img, RCoords(1,:), RCoords(2,:), ...
    'linear', extrapVal );
  out = reshape( interped, sImg );
end
