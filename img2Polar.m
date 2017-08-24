
function [outImg,rs,thetas] = img2Polar( img, varargin )
  % [outImg,rs,thetas] = img2Polar( img, [ dr, dTheta, 'method', method, ...
  %   'extrapval', extrapval ] )
  %
  % Inputs:
  % img - a 2D array
  % 
  % Optional Inputs:
  % dr - the difference between adjacent radius values (in pixels0
  % dTheta - the difference between adjacent angles (in radians)
  % method - the method of interpolation used; any method accepted by
  %   interp2 is valid
  % extrapval - the value of a pixel to output when extrapolating
  %
  % Outputs:
  % outImg - the img converted to polar coordinates with radius along the
  %   vertical dimension and angle along the horizontal dimension
  % rs - the radiuses used in polar coordinates (in pixels)
  % thetas - the angles used in polar coordinates (in radians)
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addOptional( 'dr', 1 );
  p.addOptional( 'dTheta', pi/180 );
  p.addOptional( 'maxRadius', 0 );
  p.addParameter( 'method', 'linear' );
  p.addParameter( 'extrapval', 0 );
  p.parse( varargin{:} );
  dr = p.Results.dr;
  dTheta = p.Results.dTheta;
  maxRadius = p.Results.maxRadius;
  method = p.Results.method;
  extrapval = p.Results.extrapval;


  sImg = size( img );
  if maxRadius==0, maxRadius = floor(min(size(img))/2); end;

  rs = 0:dr:maxRadius;
  thetas = -pi:dTheta:pi-dTheta;

  [rsImg,thetasImg] = ndgrid( rs, thetas );

  polarXs = rsImg .* cos(thetasImg);
  polarYs = rsImg .* sin(thetasImg);

  tmp = size2imgCoordinates( sImg );
  ys=tmp{1};  xs=tmp{2};
  [ys, xs] = ndgrid( ys, xs );

  outImg = interp2( xs, ys, img, polarXs, polarYs, method, extrapval );
end
