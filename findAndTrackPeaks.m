

function [pts1,pts2] = findAndTrackPeaks( img1, img2, varargin )
  % [pts1,pts2] = findAndTrackPeaks( img1, img2 ...
  %   [, 'N', N, 'buffer', buffer, 'w', w, 'k', k, ...
  %   'searchWidth', searchWidth, 'offset', offset ] )
  %
  % Inputs:
  % img1/img2 - 2D arrays; find features in img1 and track into img2
  % N - the maximum number of features to find in img1 (default is 50)
  % buffer - the spacing between features in img1 (default is 20)
  % featureWidth - the size of the feature
  % w - the feature width (default is 11)
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
  defaultSearch = ceil(0.2*sImg);
  defaultOffset = [0 0];
  p = inputParser;
  p.addParameter( 'N', defaultN );
  p.addParameter( 'buffer', defaultBuffer );
  p.addParameter( 'w', defaultW );
  p.addParameter( 'searchWidth', defaultSearch );
  p.addParameter( 'offset', defaultOffset );
  p.parse( varargin{:} );
  N = p.Results.N;
  buffer = p.Results.buffer;
  w = p.Results.w;
  search = p.Results.searchWidth;
  offset = p.Results.offset;

  gSig = 7;
  filtImg1 = imgaussfilt( img1, gSig );
  pts1 = findPeaks( filtImg1, N, 'buffer', buffer );
  %figure; imshowscale(img1,3); labelImgPts( 3*pts1 );  title('Image 1');

  filtImg1 = imgaussfilt( img1, 3 );
  filtImg2 = imgaussfilt( img2, 3 );
  pts2 = trackFeatures( pts1, filtImg1, filtImg2, 'kernelWidth', w, ...
    'searchWidth', search, 'offset', offset );
  %figure; imshowscale(img2,3); labelImgPts( 3*pts2 );  title('Image 2');

  maxPts2 = max( pts2, [], 2 );
  goodIndxs = find( maxPts2 > 0 );
  pts1 = pts1(goodIndxs,:);
  pts2 = pts2(goodIndxs,:);
end


