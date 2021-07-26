
function status = withinPolygon( polygon, p )
  % status = withinPolygon( polygon, p )
  %   determines if a given number of points is within a polygon
  %
  % Inputs:
  % polygon - 2D array with 2 rows and n columns consisting of polygon
  %   vertex coordinates
  % p - 2D array with 2 rows and n columns of points to be tested whether 
  % inside the polygon
  % 
  % Outputs: 
  % status - 1D array of 1s and 0s
  %   1 if point inside of polygon
  %   0 if point outside of polygon
  %
  % Written by Jolie Wang - Copyright 2021
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.
  
  q = max(polygon, [], 2); % point outside of polygon
  q = q + 1;
  nP = size( p, 2 );
  nCrosses = zeros( 1, nP );
  
  for i = 1:length(polygon)
    if i == length(polygon)
      side = [polygon(:,end) polygon(:,1)];
      sideParams = lineIntersectionParams( side, q, p );
      nCrosses = nCrosses + ( ( max( sideParams, [], 1 ) <= 1 ) & ( min( sideParams, [], 1 ) >= 0 ) );
    else
      side = [polygon(:,i) polygon(:,i+1)];
      sideParams = lineIntersectionParams( side, q, p );
      nCrosses = nCrosses + ( ( max( sideParams, [], 1 ) <= 1 ) & ( min( sideParams, [], 1 ) >= 0 ) );
    end
  end

  status = mod( nCrosses, 2 );

end

function params = lineIntersectionParams( L1, q, p )
  % L1 is a Dx2 array that specifies the line in D dimensions
  % q is a point outside of the polygon
  % p is a list of points that we would like to query

  D = size( L1, 1 );
  nP = size( p, 2 );

  pq = bsxfun( @plus, -p, q );
  dL1 = L1(:,2) - L1(:,1);  

  params = zeros( D, nP );
  b = q - L1(:,1);
  for i = 1 : nP
    A = [ dL1 pq(:,i)  ];
    params(:,i) = A \ b;
  end
end