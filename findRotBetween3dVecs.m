
function R = findRotBetween3dVecs( v1, v2 )
  % R = findRotBetweenVecs( v1, v2 )
  %
  % Inputs:
  % v1 - a 3 element array representing first vector
  % v2 - a 3 element array representing second vector
  %
  % Outputs:
  % R - a 3x3 array representing the rotation matrix
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.
  
  u1 = v1 / norm(v1,2);
  u2 = v2 / norm(v2,2);

  cp = cross( u1, u2 );
  cpM = makeCrossProdMatrix( cp );

  s = norm( cp, 2 );
  c = dot(u1,u2);  % cosine of angle between vectors
  
  R = eye(3) + cpM + cpM*cpM * (1-c)/(s*s);
end

