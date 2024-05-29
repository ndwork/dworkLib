
function [ movedPts, count ] = movePointsToSet( pts, setPts )
  % [ movedPts, samples ] = movePointsToGrid( pts, set )
  %
  % Moves points in an Real^N space to the points specified by set
  %
  % Inputs:
  % pts - a 2D array where each column corresponds to a point in R^N
  % set - a 2D array where each column corresponds to a point in R^N
  %
  % Outputs:
  % movedPts - a 2D array where each column corresponds to a moved point
  %
  % Optional Outputs:
  % count - a 1D array indicating the number of times that a pt was moved to the
  %   corresponding set point.
  %
  % Written by Nicholas Dwork, Copyright 2024
  %
  % This software is offered under the GNU General Public License 3.0.  It is offered without any warranty, 
  % expressed or implied, including the implied warranties of merchantability or fitness for a particular purpose.

  movedPts = zeros( size( pts ) );

  nSetPts = size( setPts, 2 );
  count = zeros( nSetPts, 1 );

  movedIndx = 1;
  for i = 1 : size( pts, 1 )
    distPt2Set = LpNorms( bsxfun( @minus, setPts, pts(i,:) ), 2, 1 );
    [ ~, minDistIndx ] = min( distPt2Set );
    if count( minDistIndx ) == 0
      movedPts( movedIndx, : ) = setPts( :, minDistIndx );
      movedIndx = movedIndx + 1;
    end
    count( minDistIndx ) = count( minDistIndx ) + 1;
  end

end

