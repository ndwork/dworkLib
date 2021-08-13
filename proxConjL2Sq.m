
function out = proxConjL2Sq( x, sigma, c, b )
  % out = proxConjL2Sq( x, sigma, c, b );
  %
  % Let g(x) = (c/2) || x - b ||_2^2.  This function returns the
  % proximal operator of sigma times the conjugate function of g.
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp( 'Usage:   out = proxConjL2Sq( x, sigma, c, b );' );
    if nargout > 0, out=[]; end
    return;
  end

  if nargin < 2, sigma=1; end
  if nargin < 3, c=1; end
  if nargin < 4, b=0; end

  out = ( c / ( sigma + c ) ) * ( x - sigma * b ) ;

end
