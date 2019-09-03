
function dists = findDistsBetweenPtsAndLine( pts, line )
  % dist = findDistBetweenPtsAndLine( pts, line )
  %
  % Inputs:
  % pts - a MxN real array specifying the coordinates of the points
  %   M is the dimension of the space; N is the number of points
  % line - a structure containing two elements: pt and vec
  %   line.pt is a point on the line
  %   line.vec is a vector that points in the direction of the line
  %
  % Outputs:
  % dists - a 1D array representing the shortest distance between the points and the line
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  a = line.pt;   a = a(:);
  n = line.vec;  n = n(:);

  aMinusP = bsxfun( @minus, a, pts );
  tmp = bsxfun( @times, ( aMinusP' * n )', n );
  dists = norms( aMinusP - tmp, 2, 1 );
end
