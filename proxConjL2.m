
function out = proxConjL2( x, t, b )
  % out = proxConjL2( x, t, b );
  %
  % Let g = (1/2) || x - b ||_2^2.  This function returns the proximal
  % operator of the conjugate function of t g, where t is a scalar function.
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 2, t=1; end
  if nargin < 3, b=0; end

  out = ( x - t * b ) / ( t + 1 );

end
