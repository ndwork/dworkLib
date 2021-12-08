
function out = undoRadialDistortion( img, ks, varargin )
  % out = undoRadialDistortion( img, k [, 'c', c, 'space', 'full',[] ] )
  %
  % Written according to section 7.4 of Multiple View Geometry, 2nd edition
  % by Hartley and Zisserman
  %
  % Inputs:
  % img - a 2D array or 3D array where the third dimension is color
  % ks - the radial distortion coefficients
  %   L(r) = 1 + k(1) * r + k(2) * r^2 + ...
  %
  % Optional Inputs:
  % c - the center of the image [ center_x center_y ] (default is [0; 0;])
  % space - if full, shows you all the data, otherwise, shows you area of size(img)
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
  p.addParameter( 'space', [], @(x) true );
  p.parse( varargin{:} );
  c = p.Results.c;
  space = p.Results.space;

  if numel( c ) == 0,  c = zeros(2,1); end

  sImg = size(img);
  coords = size2imgCoordinates( sImg(1:2) );

  [xs,ys] = meshgrid( coords{2}, coords{1} );

  sImg = size(img);
  if numel( space ) > 0 && strcmp( space, 'full' )
    ptUndone = applyRadialDistortion2Pts( [xs(:) ys(:)], ks, c, 'dir', -1 );
      % TODO:  Make this faster by just undoing corners
    xUndone = ptUndone(:,1);
    yUndone = ptUndone(:,2);

    minX = floor( min( xUndone(:) ) );  maxX = ceil( max( xUndone(:) ) );
    minY = floor( min( yUndone(:) ) );  maxY = ceil( max( yUndone(:) ) );
  else
    minX = 1;  maxX = sImg(2);
    minY = 1;  maxY = sImg(1);
  end

  [xOuts,yOuts] = meshgrid( minX : maxX, minY : maxY );
  sOut = [ maxY-minY+1 maxX-minX+1 size(img,3) ];

  distortedOutPts = applyRadialDistortion2Pts( [xOuts(:) yOuts(:)], ks, c );
  xDistorted = reshape( distortedOutPts(:,1), sOut(1:2) );
  yDistorted = reshape( distortedOutPts(:,2), sOut(1:2) );

  out = zeros( sOut );
  for colorIndx = 1 : size( img, 3 )
   colorChannel = img(:,:,colorIndx);
   tmp = interp2( xs, ys, colorChannel, xDistorted, yDistorted, 'linear', 0 );
   out(:,:,colorIndx) = reshape( tmp, sOut(1:2) );
  end

end



