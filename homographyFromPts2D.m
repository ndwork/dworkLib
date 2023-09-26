
function H = homographyFromPts2D( pts1, pts2 )
  % H = homographyFromPts2D( pts1, pts2 )
  %
  % Determine the homography that projects pts1 onto pts2 with the Direct
  %   Linear Transformation.
  %
  % Inputs:
  % pts1 - An Nx2 array where N is the number of points
  % pts2 - An Nx2 array where N is the number of points
  %
  % Outputs:
  % H - a 3x3 matrix representing the homography
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 2
    disp( 'Usage: H = homographyFromPts2D( pts1, pts2 )' );
    if nargout > 0, H = []; end
    return;
  end

  nPts = size( pts1, 1 );

  [newPts1,T1] = normalizePts2D( pts1 );
  [newPts2,T2] = normalizePts2D( pts2 );

  % Convert to homogeneous points
  pts1_h = ones( 3, nPts );
  pts1_h(1:2,:) = newPts1';
  pts2_h = ones( 3, nPts );
  pts2_h(1:2,:) = newPts2';

  A = zeros( 2*nPts, 9 );
  for i=1:nPts
    % First row
    A(2*i-1,4:6) = -pts2_h(3,i) * pts1_h(:,i)';
    A(2*i-1,7:9) =  pts2_h(2,i) * pts1_h(:,i)';

    % Second row
    A(2*i,1:3) =  pts2_h(3,i) * pts1_h(:,i)';
    A(2*i,7:9) = -pts2_h(1,i) * pts1_h(:,i)';
  end

  [~,~,v] = svd(A);
  h = v(:,end);
  H = reshape( h, [3 3] )';

  H = T2 \ H * T1;   % denormalize
end
