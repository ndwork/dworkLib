
function legendnice( varargin )
  % legendnice( names )
  %
  % Makes a legend with large font and puts it in the best place
  %
  % Written by Nicholas - Copyright 2017
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  legHandle = legend( varargin{:}, 'Location', 'best' );
  set( legHandle, 'FontSize', 18 )
end
