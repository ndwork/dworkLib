
function [pts1,pts2] = findAndTrackCorners( img1, img2, varargin )
  % corners = findAndTrackCorners( img1, img2 ...
  %   [, 'N', N, 'buffer', buffer, 'w', w, 'k', k, ...
  %   'searchWidth', searchWidth, 'offset', offset ] )
  %
  % Inputs:
  % img1/img2 - 2D arrays; find features in img1 and track into img2
  % N - the maximum number of features to find in img1
  % buffer - the spacing between features in img1
  % w - the feature width
  % k - Harris corner parameter (default is 0.04, set for images nominally
  %   scaled between 0 and 1)
  %
  % Optional Inputs:
  % searchWidth - the size of the search area to consider when tracking
  %   (default is 20% of the image size)
  % offset - 2 element array; the [y x] shift to center the search area
  %   (default is [0 0])
  %
  % Outputs:
  % pts1 - an Nx2 array.  The first/second column is the x/y coordinate
  %   of each feature in img1.
  % pts2 - an Nx2 array.  The first/second column is the x/y coordinate
  %   of each feature in img2.
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It 
  % is offered without any warranty expressed or implied, including the 
  % implied warranties of merchantability or fitness for a particular 
  % purpose.

  sImg = size(img1);

  defaultN = 50;
  defaultBuffer = 20;
  defaultW = 11;
  defaultK = 0.04;
  defaultSearch = ceil(0.2*sImg);
  defaultOffset = [0 0];
  p = inputParser;
  p.addParamValue( 'N', defaultN );
  p.addParamValue( 'buffer', defaultBuffer );
  p.addParamValue( 'w', defaultW );
  p.addParamValue( 'k', defaultK );
  p.addParamValue( 'searchWidth', defaultSearch );
  p.addParamValue( 'offset', defaultOffset );
  p.parse( varargin{:} );
  N = p.Results.N;
  buffer = p.Results.buffer;
  w = p.Results.w;
  k = p.Results.k;
  search = p.Results.searchWidth;
  offset = p.Results.offset;

  pts1 = findHarrisCorners(img1, N, buffer, w, k);
  %figure; imshow(img1,[]); labelImgPts( pts1 );  title('Image 1');

  pts2 = trackFeatures( pts1, img1, img2, 'kernelWidth', w, ...
    'searchWidth', search, 'offset', offset );
  %figure; imshow(img2,[]); labelImgPts( pts2 );  title('Image 2');

  maxPts2 = max( pts2, [], 2 );
  goodIndxs = find( maxPts2 > 0 );
  pts1 = pts1(goodIndxs,:);
  pts2 = pts2(goodIndxs,:);
end


