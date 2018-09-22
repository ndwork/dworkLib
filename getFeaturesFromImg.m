
function features = getFeaturesFromImg( n, varargin )
  % features = getFeaturesFromImg( n [, scale ] )
  %
  % Inputs:
  % n - the number of features to select
  % scale - the scale of the displayed image (default is 1)
  %
  % Output:
  % features - a 2 column array.  The first column are the x (or
  % horizontal) coordinates, and the second column are the y (or
  % vertical coordinates).
  %
  % Written by Nicholas Dwork - Copyright 2018
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addOptional( 'scale', 1, @isnumeric );
  p.parse( varargin{:} );
  scale = p.Results.scale;

  scaledFeatures = ginput(n);
  features = round( scaledFeatures / scale );
end