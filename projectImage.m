
function out = projectImage( img, H, varargin )
  % out = projectImageWithHomography( img, H, [ range ] )
  %
  % Inputs:
  % img - a 2D array (or a 3D array where the third dimension is color) of
  %   data to project
  % H - the 3x3 homography
  %
  % Optional Inputs:
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

  N = (range(2)-range(1)+1) * (range(4)-range(3)+1);
  xs = range(1):range(2);
  ys = range(3):range(4);
  [xs, ys] = meshgrid( xs, ys );
  coords = [ xs(:)'; ys(:)'; ones(1,N) ];

  %projCoords_h = inv(H) * coords;
  projCoords_h = H \ double( coords );
  projCoords = hom2Euc( projCoords_h );
  
  imgXs = 1:sImg(2);
  imgYs = 1:sImg(1);
  [imgXs,imgYs] = meshgrid( imgXs, imgYs );
  
  dimOut = [range(4)-range(3)+1 range(2)-range(1)+1];
  nColors = size( img, 3 );
  out = zeros( [dimOut, nColors]  );
  for i=1:nColors
    interped = interp2( imgXs, imgYs, double(img), ...
      projCoords(1,:), projCoords(2,:), 'linear', 0);
    out(:,:,i) = reshape( interped, dimOut );
  end

end
