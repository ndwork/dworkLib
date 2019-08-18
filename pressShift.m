

function pressShift
  % pressShift
  %
  % Uses the java robot to press and depress shift.
  % Called in startup so that restart exits debug mode correctly.
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  import java.awt.Robot;
  import java.awt.event.KeyEvent;

  rob = Robot;  %Create a Robot-object to do the key-pressing
  rob.keyPress( KeyEvent.VK_SHIFT );
  rob.keyRelease( KeyEvent.VK_SHIFT );
end

