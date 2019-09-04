
function out = applyRadialDistortion( img, ks, varargin )
  % out = applyRadialDistortion( img, k [, c ] )
  %
  % Written according to section 7.4 of Multiple View Geometry, 2nd edition
  % by Hartley and Zisserman
  %
  % Inputs:
  % img - a 2D array or 3D array where the third dimension is color
  % ks - the radial distortion coefficients
  %   L(r) = 1 + k(1) * r + k(2) * r^2 + ...
  % c - the center of the image [ center_x center_y ]
  %
  % Outputs:
  % out - the image with radial distortion applied
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addParameter( 'c', [], @(x) isnumeric(x) && numel(x) == 2 );
  p.parse( varargin{:} );
  c = p.Results.c;

  if numel( c ) == 0,  c = zeros(2,1); end

  sImg = size(img);
  coords = size2imgCoordinates( sImg(1:2) );

  [xs,ys] = meshgrid( coords{2}, coords{1} );
  
  distortedPts = applyRadialDistortion2Pts( [ xs(:) ys(:) ], ks, c );
  xTildes = distortedPts( :, 1 );
  yTildes = distortedPts( :, 2 );

  out = zeros( sImg );
  for colorIndx = 1 : size( img, 3 )
    colorChannel = img(:,:,colorIndx);
    tmp = interp2( xs, ys, colorChannel, xTildes, yTildes, 'linear', 0 );
    out(:,:,colorIndx) = reshape( tmp, sImg(1:2) );
  end

end
