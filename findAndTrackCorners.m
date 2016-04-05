
function [pts1,pts2] = findAndTrackCorners( img1, img2, varargin )
  % corners = findAndTrackCorners( img1, img2 ...
  %   [, 'N', N, 'buffer', buffer, 'w', w, 'k', k, 'offset', offset ] )

  defaultN = 50;
  defaultBuffer = 20;
  defaultW = 7;
  defaultK = 0.04;
  defaultOffset = [0 0];
  p = inputParser;
  p.addParamValue( 'N', defaultN );
  p.addParamValue( 'buffer', defaultBuffer );
  p.addParamValue( 'w', defaultW );
  p.addParamValue( 'k', defaultK );
  p.addParamValue( 'offset', defaultOffset );
  p.parse( varargin{:} );
  N = p.Results.N;
  buffer = p.Results.buffer;
  w = p.Results.w;
  k = p.Results.k;
  offset = p.Results.offset;

  pts1 = findHarrisCorners(img1, N, buffer, w, k);
  %showFeaturesOnImg( pts1, img1 );  title('Image 1');
  figure; imshow(img1,[]); labelImgPts( pts1 );  title('Image 1');

  pts2 = trackFeatures( pts1, img1, img2, offset );
  %showFeaturesOnImg( pts2, img2 );  title('Image 2');
  figure; imshow(img2,[]); labelImgPts( pts2 );  title('Image 2');
end

