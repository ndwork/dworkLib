
function out = scatterIntensity( x, y, z, varargin )
  % out = scatterIntensity( x, y, z, options )
  %
  % Inputs:
  % x - 1D array of horizontal values of the scatter points
  % y - 1D array of vertical values of the scatter points
  % z - 1D array of the intensities of each scatter point
  %
  % Optional Inputs:
  % options - all optional arguments of scatter
  %
  % Outputs:
  % out - handle to object
  %
  % Written by Nicholas Dwork - Copyright 2021
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  if nargin < 1
    disp( 'Usage:  out = scatterIntensity( x, y, z, varargin )' );
    return
  end
  
  if numel( x ) ~= numel( y ), error( 'x, y, and z must have the same number of elements' ); end
  if numel( x ) ~= numel( z ), error( 'x, y, and z must have the same number of elements' ); end

  pointsize = 4;
  out = scatternice( x, y, pointsize, z, varargin{:} );
end
