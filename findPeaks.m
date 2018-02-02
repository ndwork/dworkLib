
function corners = findPeaks( img, varargin )
  % corners = findPeaksCorners( img [, N, 'buffer', buffer ] )
  %
  % Inputs:
  % img - a 2D array
  % N - the number of features to identify (default is 50)
  % buffer - the minimum spacing between features (default is 20)
  %
  % Outputs:
  % corners - an Nx2 array.  The first/second column is the x/y coordinate
  %   of each feature.
  %
  % Written by Nicholas Dwork 2018
  %
  % This software is offered under the GNU General Public License 3.0.  It 
  % is offered without any warranty expressed or implied, including the 
  % implied warranties of merchantability or fitness for a particular 
  % purpose.

  defaultN = 50;
  defaultBuffer = 20;
  p = inputParser;
  p.addOptional( 'N', defaultN );
  p.addParameter( 'buffer', defaultBuffer, @isnumeric );
  p.parse( varargin{:} );
  N = p.Results.N;
  buffer = p.Results.buffer;

  sImg = size( img );

  Ix = zeros( sImg );
  Ix(:,1:end-1) = img(:,2:end) - img(:,1:end-1);
  IxSq = Ix .* Ix;

  Iy = zeros( sImg );
  Iy(1:end-1,:) = img(2:end,:) - img(1:end-1,:);
  IySq = Iy .* Iy;

  score = IxSq + IySq;

  minScore = min( score(:) );
  score(1:buffer,:) = minScore;
  score(:,1:buffer) = minScore;
  score(end-buffer:end,:) = minScore;
  score(:,end-buffer:end) = minScore;

  corners = zeros(N,2);
  for img=1:N
    [~,maxIndx] = max( score(:) );
    [y,x] = ind2sub( sImg, maxIndx );
    corners(img,1) = x;
    corners(img,2) = y;

    bIndx = max( y-buffer, 1 );
    uIndx = min( y+buffer, sImg(1) );
    LIndx = max( x-buffer, 1 );
    rIndx = min( x+buffer, sImg(2) );
    score(bIndx:uIndx,LIndx:rIndx) = minScore;
    
    if max(score(:)) == minScore, break; end;
  end
end


