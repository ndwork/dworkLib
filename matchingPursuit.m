
function x = matchingPursuit( A, b, K )
  % x = matchingPursuit( A, b, K )
  %
  % Written according to "Matching Pursuit of Images" by Bergeaud and
  % Mallat, 1995
  % performs matching pursuit optimization that greedily attempts to solve
  % minimize (1/2) || A x - b ||_2^2 subject to ||x||_0 <= k
  %
  % Note that if the set of the columns of A is orthonormal then matching
  % pursuit is the same as orthogonal matching pursuit (which generally
  % performs much better)
  %
  % Inputs:
  % A - the dictionary matrix of matching pursuit; each column is an
  %   element of the dictionary
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
  for k = 1 : K
    AHR = A' * R;  % Hermitian Transpose times b is dot product with dictionary vectors
    if k==1, x = zeros( size( AHR ) ); end
    [~,maxIndx] = max( abs( AHR ) );
    alpha = AHR( maxIndx );
    x( maxIndx ) = alpha;
    R = R - alpha * A(:,maxIndx);
    A(:,maxIndx) = 0;
  end
end
