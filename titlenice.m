
function titlenice( myTitle, varargin )
  % titlenice( myTitle, varargin )
  %
  % Inputs:
  % myTitle - String to name the graphic as well as the figure
  %
  % Optional Inputs:
  % varargin - all optional inputs accepted by title function
  %
  % Written by Nicholas Dwork - Copyright 2018
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin > 1
    title( myTitle, varargin{:} );
  else
    title( myTitle );
  end
  set( gcf, 'Name', myTitle )
  drawnow;
end
