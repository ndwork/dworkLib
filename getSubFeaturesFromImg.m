
function [features,values] = getSubFeaturesFromImg( img, n, varargin )
  % features = getSubFeaturesFromImg( img, n [, scale ] )
  %
  % Inputs:
  % n - the number of features to select
  % scale - the scale of the displayed image (default is 1)
  %
  % Output:
  % features - a 2 column array.  The first column are the x (or
  % horizontal) coordinates, and the second column are the y (or
  % vertical coordinates).
  %
  % Optional Output:
  % values - the values of the image at each feature point
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addRequired( 'img', @(x) isnumeric(x) );
  p.addRequired( 'n', @(x) isnumeric(x) && x > 0 );
  p.addOptional( 'scale', 1, @isnumeric );
  p.parse( img, n, varargin{:} );
  scale = p.Results.scale;

  mainF = figure;
  imshowscale( img, scale );

  features = zeros( n, 2 );
  subF = figure;
  for pt = 1 : n
    figure( mainF );  title([ 'Select point: ', num2str(pt) ]);
    grossPt = round( ginput( 1 ) / scale );
    hL = max( grossPt(1) - 10, 1 );
    vL = max( grossPt(2) - 10, 1 );
    hU = min( grossPt(1) + 10, size(img,2) );
    vU = min( grossPt(2) + 10, size(img,1) );
    subImg = img( vL : vU, hL : hU, : );
    figure( subF );
    subImg = subImg - min( subImg(:) );
    subImg = subImg / max( subImg(:) );
    imshowscale( subImg, scale*10 );
    features(pt,:) = double( ginput(1) ) / ( scale * 10 );
    features(pt,:) = features(pt,:) + [ hL-1 vL-1 ];
  end
  close( subF );
  close( mainF );


  if nargout > 1
    values = zeros( n, size(img,3) );
    ceiledFeatures = ceil( features );
    for fIndx = 1:n
      values(fIndx,:) = img( ceiledFeatures(fIndx,2), ceiledFeatures(fIndx,1), : );
    end
  end
end
