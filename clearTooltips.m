
function clearTooltips
  % Deletes stubborn tooltips that refuse to disappear
  %
  % Written by Nicholas - Copyright 2017
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  delete(findall(gcf,'Type','hggroup'));
end
