
function colorImg = colorImageWithMap( img, varargin )
  % colorImg = colorImageWithMap( img [, cmap ] )
  %
  % Inputs:
  % img - 2D array representing the image
  %
  % Optional Inputs:
  % cmap - the colormap (default is hot)
  %
  % Written by Nicholas - Copyright 2019
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addOptional( 'mapName', 'hot', @(x) true );
  p.parse( varargin{:} );
  mapName = p.Results.mapName;

  maxImg = max( img(:) );
  scaledImg = uint8( img / maxImg * 255 );  %#ok<NASGU>

  colorImg = 0;
  str2eval = [ 'colorImg = ind2rgb( scaledImg, ', mapName, '(256) );' ];
  eval( str2eval );
end

