
function out = falseColorize( img, varargin )
  % out = falseColorize( img [, colormapName ] )
  %
  % Inputs:
  % img - array with values between 0 and 1
  %
  % Optional Inputs:
  % colormapName - the name of the colormap to use (default is jet)
  %
  % Written by Nicholas Dwork - Copyright 2020
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addOptional( 'colormapName', 'jet', @(x) true );
  p.parse( varargin{:} );
  colormapName = p.Results.colormapName;

  eval( [ 'map = ', colormapName, '(256);' ] );

  out = ind2rgb( round(img * 255), map );
end
