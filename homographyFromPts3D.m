
function H = homographyFromPts3D( pts1, pts2 )
  % H = homographyFromPts3D( pts1, pts2 )
  %
  % Determine the homography that projects pts1 onto pts2 with the Direct
  %   Linear Transformation.
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  nPts = size( pts1, 1 );

  [newPts1,T1] = normalizePts3D( pts1 );
  [newPts2,T2] = normalizePts3D( pts2 );

  % Convert to homogeneous points
  pts1_h = ones( 4, nPts );
  pts1_h(1:3,:) = newPts1';
  pts2_h = ones( 4, nPts );
  pts2_h(1:3,:) = newPts2';

  A = zeros( 2*nPts, 16 );
  for i=1:nPts
    % First row
    A(3*i-2,1:4) = pts2_h(4,i) * pts1_h(:,i)';
    A(3*i-2,13:16) = -pts2_h(1,i) * pts1_h(:,i)';

    % Second row
    A(3*i-1,5:8) =  pts2_h(4,i) * pts1_h(:,i)';
    A(3*i-1,13:16) = -pts2_h(2,i) * pts1_h(:,i)';

    % Third row
    A(3*i,9:12) = pts2_h(4,i) * pts1_h(:,i)';
    A(3*i,13:16) = -pts2_h(3,i) * pts1_h(:,i)';
  end

  [~,~,v] = svd(A);
  h = v(:,end);
  H = reshape( h, [4 4] )';

  H = inv(T2) * H * T1;   % denormalize
end

