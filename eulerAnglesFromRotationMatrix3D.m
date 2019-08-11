
function [angles,angles2] = eulerAnglesFromRotationMatrix3D( R )
  % Finds the euler angles (phi,theta,psi) from a 3x3 rotation matrix R
  % where R(phi,theta,psi) = Rz(phi) Ry(theta) Rx(psi) according to
  % http://www.gregslabaugh.net/publications/euler.pdf
  %
  % Inputs:
  % R - a 3x3 rotation matrix
  %
  % Outputs:
  % angles - a structure containing (phi,theta,psi)
  % angles2 - unless the state is gimbal lock (theta = pi/2), then there are two
  %   valid solutions; angles2 contains the other one.
  %
  % Written by Nicholas Dwork, Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if abs( R(3,1) ) ~= 1

    theta = -asin( squeeze( R(3,1,:) ) );
    cosTheta = cos( theta );
    psi = atan2( squeeze( R(3,2,:) ) ./ cosTheta, squeeze( R(3,3,:) ) ./ cosTheta );
    phi = atan2( squeeze( R(2,1,:) ) ./ cosTheta, squeeze( R(1,1,:) ) ./ cosTheta );

    theta2 = pi - theta;
    cosTheta2 = cos( theta2 );
    psi2 = atan2( squeeze( R(3,2,:) ) ./ cosTheta2, squeeze( R(3,3,:) ) ./ cosTheta2 );
    phi2 = atan2( squeeze( R(2,1,:) ) ./ cosTheta2, squeeze( R(1,1,:) ) ./ cosTheta2 );

    angles2.phi = phi2;
    angles2.theta = theta2;
    angles2.psi = psi2;

  else

    phi = zeros(1,size(R,3));
    if R(3,1) == -1
      theta = pi/2 * ones(1,size(R,3));
      psi = atan2( squeeze( R(1,2,:) ), squeeze( R(1,3,:) ) );
    else
      theta = -pi/2;
      psi = atan2( squeeze( -R(1,2,:) ), squeeze( -R(1,3,:) ) );
    end

  end

  angles.phi = phi;
  angles.theta = theta;
  angles.psi = psi;

end
