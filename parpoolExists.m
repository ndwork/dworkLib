
function out = parpoolExists()
  % out = parpoolExists()
  %
  % Determines if a parpool has been created or not
  %
  % Written by Nicholas Dwork, Copyright 2018
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  out = numel( gcp('nocreate') );
end
