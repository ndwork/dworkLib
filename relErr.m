
function out = relErr( est, truth )
  % out = relErr( est, truth )
  %
  % Written by Nicholas Dwork - Copyright 2024
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 2
    disp( 'out = relErr( est, truth )' );
    if nargout > 0, out = []; end
    return;
  end

  out = norm( est(:) - truth(:) ) / norm( truth(:) );
end
