
function out = rotImg( img, angle, varargin )
  % out = rotImg( img, angle [, center, 'op', 'notransp'/'transp' ]  )
  %
  % Inputs:
  % img - a 2D array
  % angle - the angle to rotate in radians
  %
  % Optional Inputs:
  % center - a 2 element array specifying the center to rotate about
  % op - either 'notransp' or 'transp'
  %
  % Output:
  % out - a 2D array representing the rotated image
  %
  % Written by Nicholas Dwork - Copyright 2021
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addOptional( 'center', [0 0], @isnumeric );
  p.addParameter( 'op', 'notransp', @(x) true );
  p.parse( varargin{:} );
  center = p.Results.center;
  op = p.Results.op;

  imgCoords = size2imgCoordinates( size( img ) );
  [xs,ys] = meshgrid( imgCoords{2}, imgCoords{1} );

  pts = [ xs(:)'; ys(:)'; ];
  R = [ cos(angle) -sin(angle); sin(angle) cos(angle); ];

  Rpts = R * pts;

  out = bilinInterp2( imgCoords{2} - center(2), imgCoords{1} - center(1), ...
    img, Rpts(1,:), Rpts(2,:), 'op', op );
  out = reshape( out, size( img ) );
end
