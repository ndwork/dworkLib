
function dist = findDistBetweenPtAndLine( pt, line )
  % dist = findDistBetweenPtAndLine( pt, line )
  %
  % Inputs:
  % pt - a 1D array specifying the coordinate of the point
  % line - a structure containing two elements: pt and vec
  %   line.pt is a point on the line
  %   line.vec is a vector that points in the direction of the line
  %
  % Outputs:
  % dist - a scalar representing the shortest distance between the point and the line
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  a = line.pt;
  n = line.vec;
  aMinusP = a - pt;
  dist = norm( aMinusP - dotP( aMinusP, n ) * n );

end
