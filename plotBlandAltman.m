
function plotBlandAltman( x, y )
  % plotBlandAltman( x, y )
  %
  % Inputs:
  % x - array
  % y - array of same size as x
  %
  % Written by Nicholas Dwork - Copyright 2021
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  xyDiff = x - y;
  xyMean = reshape( mean( [ x(:) y(:) ], 2 ), size( x ) );

  figure;  scatternice( xyMean, xyDiff );
  xlabel( 'Mean', 'FontSize', 18 );
  ylabel( 'Diff', 'FontSize', 18 );

end
