
function out = makeCrossProdMatrix( v )
  % out = makeCrossProdMatrix( v )
  %
  % makes the cross product matrix out so that cross(v,u) = out * u
  %
  % Inputs:
  % v - a 3 elements array
  %
  % Outputs:
  % out - a 3x3 matrix
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.
  
  out = zeros(3);
  out(1,2) = -v(3);
  out(1,3) = v(2);
  out(2,1) = v(3);
  out(2,3) = -v(1);
  out(3,1) = -v(2);
  out(3,2) = v(1);
end
