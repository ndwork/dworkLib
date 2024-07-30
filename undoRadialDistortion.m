
function out = undoRadialDistortion( img, ks, varargin )
  % out = undoRadialDistortion( img, k [, 'c', c, 'space', 'full'/[] ] )
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

  if numel( c ) == 0,  c = [ 0; 0; ]; end

  sImg = size(img);
  nChannels = size( img, 3 );

  coords = size2imgCoordinates( sImg(1:2) );
  [xs,ys] = meshgrid( coords{2}, coords{1} );

  sImg = size(img);
  if numel( space ) > 0 && strcmp( space, 'full' )

    minY = min( ys(:) );  maxY = max( ys(:) );
    minX = min( xs(:) );  maxX = max( xs(:) );
    corners = [ [ minY, minX ]; ...
                [ minY, maxX ]; ...
                [ maxY, minX ]; ...
                [ maxY, maxX ]; ];

    ptsUndone = applyRadialDistortion2Pts( corners, ks, c, 'dir', -1 );
    %ptsRedone = applyRadialDistortion2Pts( ptsUndone, ks, c, 'dir', 1 );  % should be same as corners

    yUndone = ptsUndone(:,1);
    xUndone = ptsUndone(:,2);

    minUndoneY = floor( min( yUndone(:) ) );  maxUndoneY = ceil( max( yUndone(:) ) );
    minUndoneX = floor( min( xUndone(:) ) );  maxUndoneX = ceil( max( xUndone(:) ) );

  else
    minUndoneY = 1;  maxUndoneY = sImg(1);
    minUndoneX = 1;  maxUndoneX = sImg(2);
  end

  [ xOuts, yOuts ] = meshgrid( minUndoneX : maxUndoneX, minUndoneY : maxUndoneY );
  sOut = [ maxUndoneY - minUndoneY + 1, maxUndoneX - minUndoneX + 1, nChannels ];

  distortedOutPts = applyRadialDistortion2Pts( [ yOuts(:) xOuts(:) ], ks, c, 'dir', 1 );
  yDistorted = reshape( distortedOutPts(:,1), sOut(1:2) );
  xDistorted = reshape( distortedOutPts(:,2), sOut(1:2) );
  
  if nChannels > 1

    out = cell( 1, 1, nChannels );
    parfor colorIndx = 1 : size( img, 3 )
     colorChannel = img(:,:,colorIndx);
     tmp = interp2( xs, ys, colorChannel, xDistorted, yDistorted, 'linear', 0 );
     out{1,1,colorIndx} = reshape( tmp, sOut(1:2) );   %#ok<PFBNS>
    end
    out = cell2mat( out );

  else
    out = interp2( xs, ys, img, xDistorted, yDistorted, 'linear', 0 );
  end

end



