
function [pts1,pts2] = findAndTrackCorners( img1, img2, varargin )
  % corners = findAndTrackCorners( img1, img2 ...
  %   [, 'N', N, 'buffer', buffer, 'w', w, 'k', k, 'offset', offset ] )
  %
  % Inputs:
  % img1/img2 - 2D arrays; find features in img1 and track into img2
  % N - the number of features to find in img1
  % buffer - the spacing between features in img1
  % w - the feature width
  % k - Harris corner parameter
  % offset - 2 element array; the [y x] shift to center the search area
  %
  % Written by Nicholas Dwork 2016

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
  pts2 = trackFeatures( pts1, img1, img2, offset );
end
