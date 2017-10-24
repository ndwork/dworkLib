
function xSecsHavePassed = haveXSecsPassed( x )
  % Returns 1 if x seconds have passed since the last time this function
  %   returned 1 and returns 0 otherwise.
  % Returns 1 for the first call of the function.
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  persistent lastTime;

  if numel(lastTime) == 0
    lastTime = 0;
  end
  nowTime = cputime;

  xSecsHavePassed = ( nowTime - lastTime > x );

  if xSecsHavePassed > 0
    lastTime = nowTime;
  end
end
