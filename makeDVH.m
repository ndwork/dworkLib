
function [dvh,levels] = makeDVH( data, varargin )
  % [dvh,levels] = makeDVH( data [, 'nLevels', nLevels ] )
  %
  % splits data into equal levels and counts the number of data values
  % above each level
  %
  % Inputs:
  % data - array of values
  %
  % Optional Inputs:
  % nLevels - the number of levels
  %
  % Outputs:
  % dvh - vector of the number of data values above each level
  % levels - vector of the levels of the data expressed as the middle of 
  %   each bin
  %
  % Written by Jolie Wang - Copyright 2021
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    dvh = [];  levels = [];
    disp( 'Usage:  [dvh,levels] = makeDVH( data )' );
    return;
  end

  p = inputParser;
  p.addParameter( 'nLevels', 100, @ispositive );
  p.parse( varargin{:} );
  nLevels = p.Results.nLevels;

  dvh = zeros( nLevels, 1 );
  largest = max( data(:) );
  smallest = min( data(:) );

  dLevel = ( largest - smallest ) / nLevels;
  boundary = smallest + ( 0 : nLevels ) * dLevel;
  levels = boundary( 1 : end-1 ) + dLevel/2;

  for k = 1 : nLevels
    dvh( k ) = sum( data(:) > boundary(k) );
  end

end
