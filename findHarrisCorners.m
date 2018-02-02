

function corners = findHarrisCorners( img, varargin )
  % corners = findHarrisCorners( img [, N, 'buffer', buffer, 'w', w, 'k', k ] )
  %
  % Inputs:
  % img - a 2D array
  % N - the number of features to identify (default is 50)
  % buffer - the minimum spacing between features (default is 20)
  % w - the width of the kernel (default is 7)
  % k - Harris corner detector parameter
  %
  % Outputs:
  % corners - an Nx2 array.  The first/second column is the x/y coordinate
  %   of each feature.
  %
  % Written by Nicholas Dwork 2016
  % This software is offered under the GNU General Public License 3.0.  It 
  % is offered without any warranty expressed or implied, including the 
  % implied warranties of merchantability or fitness for a particular 
  % purpose.

  defaultN = 50;
  defaultW = 3;
  defaultBuffer = 20;
  defaultK = 0.04;
  p = inputParser;
  p.addOptional( 'N', defaultN );
  p.addParameter( 'buffer', defaultBuffer, @isnumeric );
  p.addParameter( 'w', defaultW, @isnumeric );
  p.addParameter( 'k', defaultK, @isnumeric );
  p.parse( varargin{:} );
  N = p.Results.N;
  buffer = p.Results.buffer;
  w = p.Results.w;
  k = p.Results.k;

  sImg = size( img );

  Ix = zeros( sImg );
  Ix(:,1:end-1) = img(:,2:end) - img(:,1:end-1);

  Iy = zeros( sImg );
  Iy(1:end-1,:) = img(2:end,:) - img(1:end-1,:);

  IxSq = Ix .* Ix;
  IySq = Iy .* Iy;
  IxIy = Ix .* Iy;

  %smoothFilter = fspecial( 'average', w );
  G11 = imgaussfilt( IxSq, w );
  G22 = imgaussfilt( IySq, w );
  G12 = imgaussfilt( IxIy, w );

  trG = G11 + G22;
  detG = G11 .* G22 - G12 .* G12;
  score = detG - k * trG .* trG;

  minScore = min( score(:) );
  score(1:buffer,:) = minScore;
  score(:,1:buffer) = minScore;
  score(end-buffer:end,:) = minScore;
  score(:,end-buffer:end) = minScore;

  corners = zeros(N,2);
  for i=1:N
    [~,maxIndx] = max( score(:) );
    [y,x] = ind2sub( sImg, maxIndx );
    corners(i,1) = x;
    corners(i,2) = y;

    bIndx = max( y-buffer, 1 );
    uIndx = min( y+buffer, sImg(1) );
    LIndx = max( x-buffer, 1 );
    rIndx = min( x+buffer, sImg(2) );
    score(bIndx:uIndx,LIndx:rIndx) = minScore;
    
    if max(score(:)) == minScore, break; end;
  end
end




