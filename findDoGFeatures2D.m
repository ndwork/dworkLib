
function features = findDoGFeatures2D( img, varargin )
  % features = findDoGFeatures2D( vol, ...
  %   [ 'nFeatures', nFeatures, 'dThresh', dThresh ] )

  defaultNFeatures = 100;
  defaultDThresh = 10;
  p = inputParser;
  p.addParamValue( 'nFeatures', defaultNFeatures );
  p.addParamValue( 'dThresh', defaultDThresh );
  p.parse( varargin{:} );
  nFeatures = p.Results.nFeatures;
  dThresh = p.Results.dThresh;

  fgBig = fspecial( 'gaussian', 9, 5 );
  fgSmall = fspecial( 'gaussian', 9, 3 );
  fImgBig = imfilter( img, fgBig );
  fImgSmall = imfilter( img, fgSmall );
  DoG = abs( fImgBig - fImgSmall );

  sImg = size(img);
  [y,x] = ind2sub( sImg, 1:numel(img) );
  yIndxs = reshape(y, sImg );
  xIndxs = reshape(x, sImg );

  dThreshSq = dThresh * dThresh;
  features = zeros(nFeatures,2);
  for i=1:nFeatures
    [~,maxIndx]= max( DoG(:) );
    [y,x] = ind2sub( sImg, maxIndx );
    features(i,1) = y;
    features(i,2) = x;

    xDist = abs( xIndxs - x );
    yDist = abs( yIndxs - y );
    rDistSq = xDist.*xDist + yDist.*yDist;
    DoG( rDistSq < dThreshSq ) = 0;
  end

end

