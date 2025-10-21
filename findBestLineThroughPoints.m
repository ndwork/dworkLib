
function [ pt, vec ] = findBestLineThroughPoints( pts )
  % Find the line that minimizes the L2 norm of the distances between the points and the line
  % The line is parameterized by a point and a vector
  %
  % [ pt, vec ] = findBestLineThroughPoints( pts )
  %
  % Inputs:
  % pts - an N x M array where N is the number of points and M is the dimension of the point space
  %
  % Outputs:
  % pt - an M element array specifying a point
  % vec - a normalized vector specifying the direction of the line
  %
  % Written by Nicholas Dwork - Copyright 2025
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  pt = mean( pts, 1 );  % The mean point serves as the point of the line

  centeredPts = bsxfun( @minus, pts, pt );

  [~,~,V] = svdc( centeredPts, 'econ' );
  vec = V(:,1);
end
