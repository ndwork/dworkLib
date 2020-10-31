
function out = proxConjL2Sq( x, sigma, b, c )
  % out = proxConjL2Sq( x, sigma, b, c );
  %
  % Let g(x) = (c/2) || x - b ||_2^2.  This function returns sigma times
  % the proximal operator of the conjugate function of g.
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 2, sigma=1; end
  if nargin < 3, b=0; end
  if nargin < 4, c=1; end

  out = ( x - sigma * b ) / ( sigma / c + 1 );

end
