function cam = makeOrthoCamP( scale, roll, pitch, yaw, camLoc )
  % cam = makeOrthoCamP( scale, roll, pitch, yaw, camLoc )
  %
  % Inputs:
  % scale - a scalar describing the scaling of space onto the image plane
  % roll - a positive rotation about the z axis
  % pitch - a positive rotation about the x axis
  % yaw - a positive rotation about the y axis
  % camLoc - a three element vector describing the location of the camera
  %   center (x,y,z).
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  Ry = [  cos(yaw) 0 sin(yaw); ...
                 0 1        0; ...
         -sin(yaw) 0 cos(yaw); ];

  Rx = [ 1          0          0; ...
         0 cos(pitch) -sin(pitch); ...
         0 sin(pitch)  cos(pitch); ];

  Rz = [ cos(roll) -sin(roll) 0; ...
         sin(roll)  cos(roll) 0; ...
                 0          0 1; ];

  R = Ry * Rx * Rz;

  K = [ scale 0 0; ...
        0 scale 0; ...
        0     0 1; ];

  P = [ 1 0 0 0; ...
        0 1 0 0; ...
        0 0 0 1; ];

  H = [ R -camLoc'; ...
        0 0 0 1; ];

  cam = K * P * H;

end
