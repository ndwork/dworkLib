function triangles = contoursToPolyhedron( contours, varargin )
  % triangles = contourToPolyhedron( contours, [ 'nPointsPerContour', nPointsPerContour ] )
  %
  % Inputs:
  % contours - cell array, where each element of the cell is a 2D matrix
  %   that specifies the points of the contour of the corresponding slice.
  %
  % Optional Inputs:
  % nPointsPerContour - the number of points allocated along each contour
  % 
  % Outputs:
  % triangles - a mesh of triangles that encompass the contours
  %   A 3D array of size 3 x 3 x nTriangles
  %
  % Written by Jolie Wang - Copyright 2021
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addParameter( 'nPointsPerContour', 100, @ispositive );
  p.parse( varargin{:} );
  nPointsPerContour = p.Results.nPointsPerContour;

  nContours = numel( contours );
  newContourPoints = cell( 1, nContours );
  for cIndx = 1 : nContours
    contourPts = contours{ cIndx };
    newContourPoints{ cIndx } = makeNewContourPoints( contourPts, nPointsPerContour, cIndx );
  end

  nTrianglesPerSlice = ( ( nPointsPerContour - 1 ) * 2 );
  nTriangles = ( nTrianglesPerSlice * ( nContours - 1 ) );
  triangles = zeros( 3, 3, nTriangles );

  % triangle formation pattern
  numTriangleCreated = 0;
  for j = 1 : nContours - 1
    for i = 1 : nTrianglesPerSlice
      if mod(i, 2) == 1
        point1Indx = fix(i / 2) + 1; % take from 1st contour
        point1 = newContourPoints{j}(point1Indx, :);
        point2Indx = fix(i / 2) + 1; % take from 2nd contour
        point2 = newContourPoints{j + 1}(point2Indx, :);
        point3Indx = fix(i / 2) + 2; % take from 1st contour
        point3 = newContourPoints{j}(point3Indx, :);
        point = [point1; point2; point3;];
        numTriangleCreated = numTriangleCreated + 1;
        triangles(:,:,numTriangleCreated) = point;
      else
        point1Indx = fix(i / 2); % take from 2nd contour
        point1 = newContourPoints{j + 1}(point1Indx, :);
        point2Indx = fix(i / 2) + 1; % take from 1st contour
        point2 = newContourPoints{j}(point2Indx, :);
        point3Indx = fix(i / 2) + 1; % take from 2nd contour
        point3 = newContourPoints{j + 1}(point3Indx, :);
        point = [point1; point2; point3;];
        numTriangleCreated = numTriangleCreated + 1;
        triangles(:,:,numTriangleCreated) = point;
      end
    end
  end
  
  bottomTriangles = contourToTriangles( contours{ 1 }, 1 );
  maxPoint = max(bottomTriangles);
  maxX = maxPoint(1);
  maxY = maxPoint(2);
  tempMin = norm([maxX maxY]);
  for j = 1:size(bottomTriangles, 1)
    if norm([bottomTriangles(j, 1) bottomTriangles(j, 2)]) < tempMin
      tempMin = norm([bottomTriangles(j, 1) bottomTriangles(j, 2)]);
      closestPointIndex = j;
    end
  end
  amountToShift = size(bottomTriangles, 1) - closestPointIndex;
  bottomTriangles = circshift(bottomTriangles, amountToShift, 1);
  
  topTriangles = contourToTriangles( contours{ nContours }, nContours );
  maxPoint = max(topTriangles);
  maxX = maxPoint(1);
  maxY = maxPoint(2);
  tempMin = norm([maxX maxY]);
  for j = 1:size(topTriangles, 1)
    if norm([topTriangles(j, 1) topTriangles(j, 2)]) < tempMin
      tempMin = norm([topTriangles(j, 1) topTriangles(j, 2)]);
      closestPointIndex = j;
    end
  end
  amountToShift = size(topTriangles, 1) - closestPointIndex + 1;
  topTriangles = circshift(topTriangles, amountToShift, 1);
  
  triangles = cat( 3, bottomTriangles, triangles, topTriangles );

end


function newContourPts = makeNewContourPoints( contourPts, nNewContourPoints, sliceNum )

  nIntervals = nNewContourPoints;
  nContourPts = size( contourPts, 1 );

  distBetweenContourPts = norms( circshift( contourPts, -1 ) - contourPts, 2, 2 );
  totalDistance = sum( distBetweenContourPts );
  distPerInterval = totalDistance / nIntervals;
  nIntervalsPerSegment = round( distBetweenContourPts / distPerInterval );

  newContourPts = sliceNum * ones( nNewContourPoints, 3 );
  ptIndx = 0;
  for j = 1 : nContourPts
    startSegIndx = j;
    endSegIndx = mod( startSegIndx, nContourPts ) + 1;

    distPerIntervalThisSeg = distBetweenContourPts( j ) / nIntervalsPerSegment( j );
    u = contourPts( endSegIndx, : ) - contourPts( startSegIndx, : );
    u = u / norm( u(:) ) * distPerIntervalThisSeg;

    nPtsForThisSeg = nIntervalsPerSegment(j);
    displacements = u(:) .* ( 0 : nPtsForThisSeg-1 );
    newContourPts( ptIndx + 1 : ptIndx + nPtsForThisSeg, 1:2 ) = ...
      bsxfun( @plus, displacements', contourPts( startSegIndx, : ) );

    ptIndx = ptIndx + nPtsForThisSeg;
  end

  maxs = max( newContourPts, [], 1 );
  maxX = maxs(1);  maxY = maxs(2);


  % find point with closest x value and closest y value to 0
  tempMin = norm([maxX maxY]);
  for j = 1:size(newContourPts, 1)
    if norm([newContourPts(j, 1) newContourPts(j, 2)]) < tempMin
      tempMin = norm([newContourPts(j, 1) newContourPts(j, 2)]);
      closestPointIndex = j;
    end
  end

  % reorder points
  amountToShift = size(newContourPts, 1) - closestPointIndex + 1;
  newContourPts = circshift(newContourPts, amountToShift, 1);
end


function faceTriangles = contourToTriangles( contourPoints, sliceNum )

  x = contourPoints(:, 1);
  y = contourPoints(:, 2);
  polyin = polyshape({x}, {y});
  T = triangulation( polyin );
  faceTriangles = zeros(3, 3, size(T.ConnectivityList, 1));
  for i = 1:size(faceTriangles, 3)
    point1Indx = T.ConnectivityList(i, 1, :);
    point1 = [T.Points(point1Indx, :) sliceNum];
    point2Indx = T.ConnectivityList(i, 2, :);
    point2 = [T.Points(point2Indx, :) sliceNum];
    point3Indx = T.ConnectivityList(i, 3, :);
    point3 = [T.Points(point3Indx, :) sliceNum];
    point = [point1; point2; point3;];
    faceTriangles(:,:,i) = point;
  end

end
