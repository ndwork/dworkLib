
function pts2 = trackFeatures( pts1, img1, img2, varargin )
  % pts2 = trackFeatures( pts1, img1, img2 [, 'searchWidth', searchWidth, ...
  %   'kernelWidth', kernelWidth, 'offset', offset] )
  %
  % Algorithm tracks points from img1 into img2 using Normalized Cross
  % correlation according to "Fast Normalized Cross Correlation" by Lewis.
  %
  % Inputs:
  % pts1 - an Nx2 array specifying the x/y location of each feature in img1
  %   It has size Nx2 where N is the number of features
  % img1 - a 2D array representing the first image
  % img2 - a 2D array representing the second image
  %   features will be tracked into this image
  %
  % Optional Inputs:
  % searchWidth - a 2 element array specifying the size of the search area
  %   by default, the search area is 20% of the image size
  % kernelWidth - a 2 element array specifying the size of the kernel
  % offset - a 2 element array specifying the y/x shift to center the
  %   search area (default is [0 0])
  %
  % Output:
  % pts2 - an Nx2 array specifying the x/y location of each feature in img2
  %
  % This software is offered under the GNU General Public License 3.0.  It 
  % is offered without any warranty expressed or implied, including the 
  % implied warranties of merchantability or fitness for a particular 
  % purpose.

  sImg = size( img1 );
  nPts = size( pts1, 1 );

  defaultSearch = ceil( 0.2 * sImg );
  defaultKernel = 11;
  p = inputParser;
  p.addParameter( 'searchWidth', defaultSearch );
  p.addParameter( 'kernelWidth', defaultKernel );
  p.addParameter( 'offset', [0 0] );
  p.parse( varargin{:} );
  offset = p.Results.offset;
  kw = p.Results.kernelWidth;
  sw = p.Results.searchWidth;

  if numel(kw)==1, kw=[kw kw]; end;
  if numel(sw)==1, sw=[sw sw]; end;
  hkw = floor( kw/2 );  % half kernel width
  hsw = floor( sw/2 );  % half search width

  pts2x = zeros(nPts,1);
  pts2y = zeros(nPts,1);
  parfor i=1:nPts
    x = pts1(i,1);
    y = pts1(i,2);
    %figure; imshowscale( img1 ); title( 'img1' );  labelImgPts( [x,y] );
    %figure; showFeaturesOnImg( [x,y], img1 );  title('img1');
    if y-hkw(1) < 1 || y+hkw(1) > sImg(1) || ...
       x-hkw(2) < 1 || x+hkw(2) > sImg(2), ...
      continue;
    end;
    template = img1(y-hkw(1):y+hkw(1),x-hkw(2):x+hkw(2));
    %figure; imshowscale( template, 5 ); title( 'template' );

    sB = max( y-hsw(1)+offset(1), 1 );  % search bottom
    sT = min( y+hsw(1)+offset(1), sImg(1) );  % search top
    sL = max( x-hsw(2)+offset(2), 1 );  % search left
    sR = min( x+hsw(2)+offset(2), sImg(2) );  % search right
    search = img2(sB:sT,sL:sR);
    if numel(search) == 0, continue; end;
    if min(size(search) - size(template)) < 2, continue; end;

    ncc = normxcorr2(template, search);
    sSearch = size(search);
    %figure; imshowscale( search ); title( 'search' );
    ncc = cropData( ncc, sSearch );
    [maxNCC,maxNccIndx] = max( ncc(:) );                                                   %#ok<ASGLU>
    [nccY,nccX] = ind2sub( sSearch, maxNccIndx );
    %figure; imshowscale( ncc ); title( 'ncc' );  labelImgPts( [nccX,nccY] );
    %figure; showFeaturesOnImg( [nccX,nccY], ncc );
    pts2x(i) = nccX + sL - 1;
    pts2y(i) = nccY + sB - 1;
    %figure; imshowscale( img2 ); title( 'img2' );  labelImgPts( [pts2x(i),pts2y(i)] );
    %figure; showFeaturesOnImg( [pts2x(i),pts2y(i)], img2 );
    
    %close all;
  end
  pts2 = [ pts2x pts2y ];

end
