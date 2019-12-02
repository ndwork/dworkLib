
function x = orthogonalMatchingPursuit( A, b, K )
  % x = orthogonalMatchingPursuit( A, b, K )
  %
  % Written according to "Signal Reconvery From Random Measurements Via
  % Orthogonal Matching Pursuit" by Tropp and Gilbert and "Orthogonal
  % Matching Pursuit: Recursive Function Approximation with Applications to
  % Wavelet Decomposition" by Pati and Krishnaprasad
  %
  % Written by Nicholas Dwork, Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if K > max( size( A ) )
    error( 'K must be less than the number of columns or rows in A' );
  end

  R = b;
  newA = zeros( size(A,1), K );
  indxs = zeros( K, 1 );

  for k = 1 : K
    AHR = A' * R;  % Hermitian Transpose times b is dot product with dictionary vectors
    [~,maxIndx] = max( abs( AHR ) );
    indxs( k ) = maxIndx;

    newA(:,k) = A(:,maxIndx);
    A(:,maxIndx) = 0;
    xValues = newA(:,1:k) \ b;

    bHat = newA(:,1:k) * xValues;
    R = b - bHat;
  end

  x = zeros( size(A,2), 1 );
  x( indxs ) = xValues;
end
