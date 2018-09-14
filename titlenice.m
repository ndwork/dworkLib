
function titlenice( myTitle )
  % titlenice( myTitle )
  %
  % Inputs:
  % myTitle - String to name the graphic as well as the figure
  %
  % Written by Nicholas Dwork - Copyright 2018
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  title( myTitle );
  set( gcf, 'Name', myTitle )
end
