
function out = proxConjHuber( x, sigma, t, mu, b )
  % out = proxConjHuber( x, sigma, t, mu )
  %
  % Calculates the proximal operator of sigma times the conjugate of f(x) = t * huber( x - b, mu ).
  % It uses the Moreau decomposition to calculate the proximal operator of the conjugate.
  %
  % Inputs:
  % x - an
  % t - a scalar
  % mu - the Huber penalty parameter
  % b - a scalar or array representing a translation of the input
  %
  % Written by Nicholas Dwork - Copyright 2021
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  if nargin < 1
    disp( 'Usage:  out = proxHuber( x, t, mu )' );
    if nargout > 0, out = []; end
    return
  end

  if nargin < 2, sigma = 1; end
  if nargin < 3, t = 1; end
  if nargin < 4, mu = 1; end

  if nargin < 5  % b = 0;
    out = x - sigma * proxHuber( x / sigma , t / sigma, mu );

  else
    out = x - sigma * proxHuber( x / sigma , t / sigma, mu, b );

  end

end

