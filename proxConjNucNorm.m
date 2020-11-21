
function out = proxConjNucNorm( in, sigma, t )
  % out = proxConjNucNorm( in, t )
  %
  % Returns the proximal operator of the conjugate of f(X) = t || X ||_*
  % where X is the input matrix.
  %
  % Inputs:
  % in - an input matrix
  % t - the thresholding value
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 2, sigma = 1; end
  if nargin < 3, t = 1; end

  if sigma == 0
    out = in;
    return;
  end

  [u,s,v] = svd( in, 'econ' );

  sPrime = min( diag(s), t );  % singular values are always real and non-negative
  %sPrime = proxConjL1( diag(s), sigma, t );

  out = u * diag( sPrime ) * v';
end
