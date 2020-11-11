
function out = figureExists()
  % out = figureExists()
  %
  % Returns true if a figure exists and false otherwise
  %
  % Written by Nicholas Dwork - Copyright 2020
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  g = groot;
  out = ~isempty( g.Children );

end
