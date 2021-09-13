
function plotBlandAltman( x, y, varargin )
  % plotBlandAltman( x, y [, 'range', range ] )
  %
  % Inputs:
  % x - array
  % y - array of same size as x
  %
  % Optional Inputs:
  % range - either 'nice' or an argument that is passed into axis command
  %   By default, all data is shown
  %
  % Written by Nicholas Dwork - Copyright 2021
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  p = inputParser;
  p.addParameter( 'range', [] );
  p.addParameter( 'thresh', 0.95 );
  p.parse( varargin{:} );
  range = p.Results.range;
  thresh = p.Results.thresh;

  xyMean = reshape( mean( [ x(:) y(:) ], 2 ), size( x ) );

  xyDiff = x - y;
  meanDiff = mean( xyDiff(:) );
  stdDiff = std( xyDiff(:) );

  scatternice( xyMean, xyDiff );
  if numel( range ) ~= 0
    if strcmp( range, 'nice' )
      xyMeanLow = findValueBelowFraction( xyMean(:), thresh );
      xyMeanHigh = findValueBelowFraction( xyMean(:), 1-thresh );

      xyDiffLow = findValueBelowFraction( xyDiff(:), thresh );
      xyDiffHigh = findValueBelowFraction( xyDiff(:), 1-thresh );

      xyDiffLow = min( xyDiffLow, meanDiff - 1.2 * stdDiff );      
      xyDiffHigh = max( xyDiffHigh, meanDiff + 1.2 * stdDiff );

      axis( [ xyMeanLow xyMeanHigh xyDiffLow xyDiffHigh ] );
    else
      axis( range );
    end
  end

  xl = xlim;
  hold on;
  plotnice( xl, [ meanDiff, meanDiff ], 'b' );
  plotnice( xl, [ meanDiff + stdDiff, meanDiff + stdDiff ], 'r--' );
  plotnice( xl, [ meanDiff - stdDiff, meanDiff - stdDiff ], 'r--' );

  xlabel( 'Mean', 'FontSize', 18 );
  ylabel( 'Diff', 'FontSize', 18 );

end
