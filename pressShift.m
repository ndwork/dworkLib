

function pressShift
  import java.awt.Robot;
  import java.awt.event.KeyEvent;

  rob = Robot;  %Create a Robot-object to do the key-pressing
  rob.keyPress( KeyEvent.VK_SHIFT );
  rob.keyRelease( KeyEvent.VK_SHIFT );
end

