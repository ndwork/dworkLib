
function out = rect( x, a )
  % out = rect( x, a )
  %
  % outputs rect( x / a );
  %
  % Inputs:
  % x - array of unscaled domain values
  % a - either a scalar or array same size as x
  %
  % Outputs:
  % out - array of rect values the same size as x
  %
  % Written by Nicholas Dwork - Copyright 2021
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp( 'Usage:  out = rect( x, a )' );
    if nargout > 0, out = []; end
    return;
  end

  out = abs( x ./ a ) < 0.5;
end
