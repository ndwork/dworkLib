
function [H,inlierIndxs,outlierIndxs] = ransacDltHomographyFromPts2D( ...
  pts1, pts2, thresh )
  % [H,inlierIndxs,outlierIndxs] = ransacDltHomographyFromPts2D( ...
  %   pts1, pts2, thresh )
  % This function finds the 3x3 homography H such that pt2 = H * pt1;
  % It implements RANSAC with the Direct Linear Transformation algorithm.
  %
  % Inputs:
  % pts1/pts2 - 2D arrays of size Mx2
  %    M is the number of points
  % thresh - points within this distance (in pixels) are considered inliers
  %
  % Outputs:
  % H - the 3x3 Homography
  %
  % Algorithm 4.5 of Multiple View Geometry, 2nd ed by Hartley and
  %   Zisserman
  % Written by Nicholas Dwork 2016

  nPts = size( pts1, 1 );
  p = 0.99;

  % Convert to homogeneous coordinates
  pts1_h = ones( 3, nPts );
  pts1_h(1:2,:) = pts1';

  sample_count = 0;
  bestNInliers = 0;
  nRansac = 1500;
  while nRansac > sample_count
    disp(['Working on ', num2str(sample_count), ' of ', num2str(nRansac)]);
    indxs = randsample(nPts,4);
    subset1 = pts1(indxs,:);
    subset2 = pts2(indxs,:);
    H = dltHomographyFromPts2D( subset1, subset2 );

    % determine the number of inliers
    proj1_h = H * pts1_h;
    proj1 = hom2Euc( proj1_h )';

    diffs = pts2 - proj1;
    dists = sqrt( diffs(:,1).*diffs(:,1) + diffs(:,2).*diffs(:,2) );

    nInliers = sum( dists < thresh );
    if nInliers > bestNInliers
      bestNInliers = nInliers;
      inlierIndxs = find( dists <= thresh );
      outlierIndxs = find( dists > thresh );

      epsilon = 1 - nInliers / nPts;
      nRansac = log( 1 - p ) / log( 1 - (1-epsilon)^4 );
    end

    sample_count = sample_count + 1;
  end

  H = dltHomographyFromPts2D( pts1(inlierIndxs,:), pts2(inlierIndxs,:) );
end
