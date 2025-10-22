
function out = applyRadialDistortion( img, ks, varargin )
  % out = applyRadialDistortion( img, k [, 'c', c ] )
  %
  % Written according to section 7.4 of Multiple View Geometry, 2nd edition
  % by Hartley and Zisserman
  %
  % Inputs:
  % img - a 2D array or 3D array where the third dimension is color
  % ks - the radial distortion coefficients
  %   L(r) = 1 + k(1) * r + k(2) * r^2 + ...
  % c - the origin of the distortion [ center_x center_y ]
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
  p.addParameter( 'dir', [], @(x) numel(x) == 0 || x == 1 || x == -1 );
  p.addParameter( 'space', [], @(x) true );
  p.parse( varargin{:} );
  c = p.Results.c;
  dir = p.Results.dir;
  space = p.Results.space;

  if numel( c ) == 0,  c = zeros(2,1); end

  sImg = size( img );
  coords = size2imgCoordinates( sImg(1:2) );
  [xs,ys] = meshgrid( coords{2}, coords{1} );

  sOut = sImg;
  if numel( space ) > 0 && strcmp( space, 'full' )

    ptsHat = applyRadialDistortion2Pts( [ xs(:) ys(:) ], ks, c, 'dir', -dir );
    xHats = ptsHat(:,1);
    yHats = ptsHat(:,2);

    minHatY = floor( min( yHats(:) ) );   maxHatY = ceil( max( yHats(:) ) );
    minHatX = floor( min( xHats(:) ) );   maxHatX = ceil( max( xHats(:) ) );
    [ xDestination, yDestination ] = meshgrid( minHatX : maxHatX, minHatY : maxHatY );

    srcPts = applyRadialDistortion2Pts( [ xDestination(:) yDestination(:) ], ks, c, 'dir', dir );
    xOuts = reshape( srcPts(:,1), size(xDestination) );
    yOuts = reshape( srcPts(:,2), size(yDestination) );

    sOut = sImg;
    sOut(1:2) = size( xOuts );

  else
    ptsTilde = applyRadialDistortion2Pts( [ xs(:) ys(:) ], ks, c, 'dir', dir );
    xTildes = reshape( ptsTilde( :, 1 ), sImg(1:2) );
    yTildes = reshape( ptsTilde( :, 2 ), sImg(1:2) );

    xOuts = xTildes;
    yOuts = yTildes;
  end

  out = zeros( sOut );
  for colorIndx = 1 : size( img, 3 )
    colorChannel = img(:,:,colorIndx);
    tmp = interp2( xs, ys, colorChannel, xOuts, yOuts, 'linear', 0 );
    out(:,:,colorIndx) = tmp;
  end

end
