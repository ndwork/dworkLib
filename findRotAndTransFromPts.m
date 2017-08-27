
function [R,t] = findRotAndTransFromPts( pts1, pts2 )
  % [R,t] = findRotAndTrans( pts1, pts2 )
  % This function finds the rotation matrix R and the translation matrix t
  %   such that pt2 = R * pt1 + t;
  % It implements the Kabsch algorithm
  % https://en.wikipedia.org/wiki/Kabsch_algorithm
  %
  % Inputs:
  % pts1/pts2 - 2D arrays of size MxN
  %    M is the number of points
  %    N is the number of dimensions
  %
  % Outputs:
  % R - the NxN rotation matrix
  % t - the N element translation vector
  %
  % Written by Nicholas Dwork

  [M,N] = size( pts1 );

  mean1 = mean( pts1, 1 );
  mean2 = mean( pts2, 1 );
  shift1 = pts1 - repmat( mean1, [M 1] );
  shift2 = pts2 - repmat( mean2, [M 1] );

  A = shift1' * shift2;
  [u,~,v] = svd(A);

  d = sign( det( v*u' ) );

  R = v * diag([ones(N-1,1) d]) * u';
  t = mean2' - R*mean1';
end
