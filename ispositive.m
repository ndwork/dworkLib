
function out = ispositive( x )
  % out = ispositive( x )
  %
  % Inputs:
  % x - a scalar or array
  %
  % Output:
  % true if all elements of x are positive and false otherwise
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  out = isnumeric(x) && min( x(:) > 0 );
end
