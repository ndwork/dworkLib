
function out = haveXSecsPassed( x )
  % Returns 1 if x seconds have passed; returns 0 otherwise
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  persistent thisTime;

  if numel(thisTime) == 0
    lastTime = 0;
  else
    lastTime = thisTime;
  end
  thisTime = tic;
  
  out = thisTime - lastTime > x;
end
