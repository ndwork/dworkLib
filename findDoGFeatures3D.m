
function features = findDoGFeatures3D( vol, varargin )
  % features = findDifferenceOfGaussianFeatures3D( vol, ...
  %   [ 'nFeatures', nFeatures, 'dThresh', dThresh ] )
  % nFeatures is the desired number of features.  May return fewer than
  %   nFeatures if image is to small.
  % dThresh is the minimum distance two features can be within each other

  defaultNFeatures = 100;
  defaultDThresh = 10;
  p = inputParser;
  p.addParamValue( 'nFeatures', defaultNFeatures );
  p.addParamValue( 'dThresh', defaultDThresh );
  p.parse( varargin{:} );
  nFeatures = p.Results.nFeatures;
  dThresh = p.Results.dThresh;

  fBig = makeGaussFilter3D( 9, 5 );
  fSmall = makeGaussFilter3D( 9, 3 );
  gVol5 = imfilter(vol, fBig);
  gVol3 = imfilter(vol, fSmall);
  DoG = abs( gVol5 - gVol3 );

  sVol = size(vol);
  [y,x,z] = ind2sub( sVol, 1:numel(vol) );
  yIndxs = reshape(y, sVol );
  xIndxs = reshape(x, sVol );
  zIndxs = reshape(z, sVol );

  dThreshSq = dThresh * dThresh;
  features = zeros(nFeatures,3);
  for i=1:nFeatures
    [~,maxIndx]= max( DoG(:) );
    [y,x,z] = ind2sub( sVol, maxIndx );
    features(i,1) = y;
    features(i,2) = x;
    features(i,3) = z;
    
    xDist = abs( xIndxs - x );
    yDist = abs( yIndxs - y );
    zDist = abs( zIndxs - z );
    rDistSq = xDist.*xDist + yDist.*yDist + zDist.*zDist;
    DoG( rDistSq < dThreshSq ) = 0;
    
    if max( DoG ) == 0, break; end;
  end

  features = features(1:i,:);
end

