
function H = dltHomographyFromPts2D( pts1, pts2 )
  % H = dltHomographyFromPts2D( pts1, pts2 )
  %
  % Determine the homography that projects pts1 onto pts2 with the Direct
  %   Linear Transformation.
  %
  % Written by Nicholas Dwork 2016

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

  H = inv(T2) * H * T1;   % denormalize
end
