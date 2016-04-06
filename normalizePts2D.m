
function [newPts,T] = normalizePts2D( pts )
  % Inputs:
  % pts - 2D array of size Nx2 where N is the number of points
  %   First/second column is x/y locations
  %
  % Outputs:
  % newPts - pts scaled and translated so that centroid is at 0 and average
  %   distance from center is sqrt(2)
  %
  % Written by Nicholas Dwork 2016

  t = -mean( pts, 1 );
  
  tPts = zeros( size(pts) );
  tPts(:,1) = pts(:,1) + t(1);
  tPts(:,2) = pts(:,2) + t(2);

  tDists = sqrt( tPts(:,1).*tPts(:,1) + tPts(:,2).*tPts(:,2) );
  meanDist = mean( tDists );

  scale = sqrt(2) / meanDist;
  newPts = tPts * scale;

  if nargout > 1
    T1 = eye(3);
    T1(1:2,3) = t;
    T2 = diag([scale scale 1]);
    T = T2 * T1;
  end
  
end
