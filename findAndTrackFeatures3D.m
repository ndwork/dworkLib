
function [pts1, pts2] = findAndTrackFeatures3D( vol1, vol2, varargin )
  % corners = findAndTrackCorners( vol1, vol2 ...
  %   [, 'N', N, 'buffer', buffer, 'w', w, ...
  %   'searchWidth', searchWidth, 'offset', offset ] )
  %
  % Inputs:
  % vol1/vol2 - 2D arrays; find features in img1 and track into img2
  % N - the maximum number of features to find in img1
  % buffer - the spacing between features in img1
  % w - the feature width
  %
  % Optional Inputs:
  % searchWidth - the size of the search area to consider when tracking
  %   (default is 20% of the image size)
  % offset - 2 element array; the [y x z] shift to center the search area
  %   (default is [0 0 0])
  %
  % Written by Nicholas Dwork - Copyright 2016


  sVol = size(vol1);

  defaultN = 50;
  defaultBuffer = 20;
  defaultW = 11;
  defaultSearch = ceil(0.2*sVol);
  defaultOffset = [0 0 0];
  p = inputParser;
  p.addParamValue( 'N', defaultN );
  p.addParamValue( 'buffer', defaultBuffer );
  p.addParamValue( 'w', defaultW );
  p.addParamValue( 'searchWidth', defaultSearch );
  p.addParamValue( 'offset', defaultOffset );
  p.parse( varargin{:} );
  N = p.Results.N;
  buffer = p.Results.buffer;
  w = p.Results.w;
  search = p.Results.searchWidth;
  offset = p.Results.offset;


  pts1 = findDoGFeatures3D( vol1, 'nFeatures', N, 'buffer', buffer );

  pts2 = trackFeatures3D( pts1, vol1, vol2, 'kernelWidth', w, ...
    'searchWidth', search, 'offset', offset );

  maxPts2 = max( pts2, [], 2 );
  goodIndxs = find( maxPts2 > 0 );
  pts1 = pts1(goodIndxs,:);
  pts2 = pts2(goodIndxs,:);

end
