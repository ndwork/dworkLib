
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
  
  numPoints = length(polygon);
  slice = [numPoints numPoints - 1];
  
  for i = 3:numPoints
    if mod(i, 2) == 1
      slice(i) = numPoints - slice(i - 1);
    else
      slice(i) = (numPoints - 1) - slice(i - 1);
    end
  end
  
  fullStatus = [];
  for i = 1:numPoints - 2
    triangle = [polygon(:,slice(i)) polygon(:,slice(i+1)) polygon(:,slice(i+2))];
    oneStatus = withinTriangle( triangle, p );
    fullStatus = [fullStatus; oneStatus]; %#ok<AGROW>
  end
  
  status = sum(fullStatus);

end
