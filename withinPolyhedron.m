
function status = withinPolyhedron( triangles, p )
  % status = withinPolyhedron( polyhedron, p )
  %   determines if a given number of points is within a polyhedron
  %
  % Inputs:
  % triangles - concatenated n triangles that each have 3 rows and 3 
  %   columns, which are polyhedron vertex coordinates of triangular facets
  % p - array with n rows and 3 columns of points to be tested whether 
  % inside the polyhedron
  % 
  % Outputs: 
  % status - 1D array of 1s and 0s
  %   1 if point inside of polyhedron
  %   0 if point outside of polyhedron
  %
  % Written by Jolie Wang - Copyright 2021
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  tiny = 1d-7;
  dimensions = size( triangles );
  numTriangles = dimensions( 3 );
  numPoints = size( p, 1 );
  numCrosses = zeros( 1, numPoints );
  eachMax = max(triangles) + 1;
  xValue = max([eachMax(:,1,1) eachMax(:,1,2) eachMax(:,1,3) eachMax(:,1,4)]);
  yValue = max([eachMax(:,2,1) eachMax(:,2,2) eachMax(:,2,3) eachMax(:,2,4)]);
  zValue = max([eachMax(:,3,1) eachMax(:,3,2) eachMax(:,3,3) eachMax(:,3,4)]);
  q = [xValue yValue zValue]; % point outside of polyhedron
  
  for j = 1:numPoints
    for i = 1:numTriangles
      a = triangles(1,:,i);
      b = triangles(2,:,i);
      c = triangles(3,:,i);
      
      % normal vector
      n = cross(b - a, c - a);
      if abs( dotP( n, p(j,:) - a ) ) < tiny
        % point lies on the boundary
        numCrosses(j) = 0;
        break; 
      end
      [ phi, theta1, ~ ] = cart2sph( n(1), n(2), n(3) );
      theta = pi/2 - theta1;

      % translate a to the origin
      b = b - a;
      c = c - a;
      tempQ = q - a;
      tempP = p(j,:) - a;
      a = a - a;

      % rotate so triangle lies in xy plane
      b = b * rotz( rad2deg(phi) ) * roty( rad2deg(theta) );
      c = c * rotz( rad2deg(phi) ) * roty( rad2deg(theta) );
      tempQ = tempQ * rotz( rad2deg(phi) ) * roty( rad2deg(theta) );
      tempP = tempP * rotz( rad2deg(phi) ) * roty( rad2deg(theta) );

      % line segment intersect triangle
      if (tempP(3) > 0 && tempQ(3) < 0) || (tempP(3) < 0 && tempQ(3) > 0)
        if tempP(3) ~= tempQ(3)
          t = -tempP(3) / (tempQ(3) - tempP(3));
          x = tempP(1) + (tempQ(1) - tempP(1)) * t;
          y = tempP(2) + (tempQ(2) - tempP(2)) * t;
          xy = [ x y ]';
          shape = [ [a(1) b(1) c(1)]; [a(2) b(2) c(2)]; ];
          in = withinPolygon( shape, xy );
          if in == 1
            numCrosses(j) = numCrosses(j) + 1;
          end
        end
      end
    end
  end

  status = mod(numCrosses, 2);
end
