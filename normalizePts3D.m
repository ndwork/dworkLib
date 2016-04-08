
function [newPts,T] = normalizePts3D( pts )
  % Inputs:
  % pts - 2D array of size Nx3 where N is the number of points
  %   First/second/third column is x/y/z locations
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
  tPts(:,3) = pts(:,3) + t(3);

  tDists = sqrt( tPts(:,1).*tPts(:,1) + tPts(:,2).*tPts(:,2) + ...
    tPts(:,3).*tPts(:,3) );
  meanDist = mean( tDists );

  scale = sqrt(2) / meanDist;
  newPts = tPts * scale;

  if nargout > 1
    T1 = eye(4);
    T1(1:3,4) = t;
    T2 = diag([scale scale scale 1]);
    T = T2 * T1;
  end

end
