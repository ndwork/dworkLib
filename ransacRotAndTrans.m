
function [R,t,inlierIndxs] = ransacRotAndTrans( pts1, pts2, thresh )
  % [R,t,inlierIndxs] = ransacRotAndTrans( pts1, pts2, thresh )
  % This function finds the rotation matrix R and the translation matrix t
  %   such that pt2 = R * pt1 + t;
  % It implements RANSAC with the Kabsch algorithm
  % https://en.wikipedia.org/wiki/Kabsch_algorithm
  %
  % Inputs:
  % pts1/pts2 - 2D arrays of size MxN
  %    M is the number of points
  %    N is the number of dimensions
  % thresh - points within this distance are considered inliers
  %
  % Outputs:
  % R - the NxN rotation matrix
  % t - the N element translation vector
  % inlierIndxs - vector indicating which points were used in final
  %   determination of transformation
  %
  % Algorithm 4.5 of Multiple View Geometry, 2nd ed by Hartley and
  %   Zisserman
  % Written by Nicholas Dwork
  %
  % This software is offered under the GNU General Public License 3.0.  It 
  % is offered without any warranty expressed or implied, including the 
  % implied warranties of merchantability or fitness for a particular 
  % purpose.

  [M, N] = size( pts1 );
  p = 0.99;

  sample_count = 0;
  bestNInliers = 0;
  nRansac = 1500;
  while nRansac > sample_count
    %disp(['Working on ', num2str(sample_count), ' of ', num2str(nRansac)]);
    indxs = randsample(M,N);  % N points determine R,t
    subset1 = pts1(indxs,:);
    subset2 = pts2(indxs,:);
    [R,t] = findRotAndTrans( subset2, subset1 );

    % determine the number of inliers
    rot2 = R * transpose(pts2);
    aligned2 = rot2';
    aligned2 = aligned2 + repmat( transpose(t), [M 1] );

    diffs = pts1 - aligned2;
    dists = sqrt( diffs(:,1).*diffs(:,1) + diffs(:,2).*diffs(:,2) );

    nInliers = sum( dists < thresh );
    if nInliers > bestNInliers
      bestNInliers = nInliers;
      inlierIndxs = find( dists < thresh );

      epsilon = 1 - nInliers / M;
      nRansac = log( 1 - p ) / log( 1 - (1-epsilon)^N );
    end

    sample_count = sample_count + 1;
  end

  [R,t] = findRotAndTrans( pts1(inlierIndxs,:), pts2(inlierIndxs,:) );
end
