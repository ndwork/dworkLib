
function status = withinTriangle( triangle, p )
  % status = withinTriangle( triangle, p )
  %   determines if a given number of points is within a triangle
  %
  % Inputs:
  % triangle - 2D array with 2 rows and n columns consisting of triangle
  % vertex coordinates
  % p - 2D array with 2 rows and n columns of points to be tested whether 
  % inside the triangle
  % 
  % Outputs: 
  % status - 1D array of 1s and 0s
  %   1 if point inside of triangle
  %   0 if point outside of triangle
  %
  % Written by Nicholas Dwork - Copyright 2021
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.
  
  a = triangle(:,1);
  b = triangle(:,2);
  c = triangle(:,3);
  q = b + (b - a);  % point outside of the triangle

  ac = [ a c ];   acParams = lineIntersectionParams( ac, q, p );
  bc = [ b c ];   bcParams = lineIntersectionParams( bc, q, p );

  nP = size( p, 2 );
  nCrosses = zeros( 1, nP );
  nCrosses = nCrosses + ( ( max( acParams, [], 1 ) <= 1 ) & ( min( acParams, [], 1 ) >= 0 ) );
  nCrosses = nCrosses + ( ( max( bcParams, [], 1 ) <= 1 ) & ( min( bcParams, [], 1 ) >= 0 ) );

  status = mod( nCrosses, 2 );
end


function params = lineIntersectionParams( L1, q, p )
  % L1 is a Dx2 array that specifies the line in D dimensions
  % q is a point outside of the triangle
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

