
function [outImg,rs,thetas] = img2Polar( img, varargin )
  % [outImg,rs,thetas] = img2Polar( img, [ 'dr', dr, 'dTheta', dTheta, ...
  %   'rs', rs, 'thetas', thetas, 'method', method, 'extrapval', extrapval ] )
  %
  % Inputs:
  % img - a 2D array
  % 
  % Optional Inputs:
  % dr - the difference between adjacent radius values (in pixels)
  % dTheta - the difference between adjacent angles (in radians)
  % rs - a vector specifying the radius values (in pixels)
  %      if specified, dr is ignored
  % thetas - a vector specifying the angles (in radians)
  %      if specified, dTheta is ignored
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

  defaultDr = 1;
  defaultDTheta = pi/180;
  defaultRs = [];
  defaultThetas = [];
  p = inputParser;
  p.addOptional( 'dr', defaultDr, @isnumeric );
  p.addOptional( 'dTheta', defaultDTheta, @isnumeric );
  p.addOptional( 'maxRadius', 0, @isnumeric );
  p.addParameter( 'rs', defaultRs );
  p.addParameter( 'thetas', defaultThetas );
  p.addParameter( 'method', 'linear' );
  p.addParameter( 'extrapval', 0, @isnumeric );
  p.parse( varargin{:} );
  dr = p.Results.dr;
  dTheta = p.Results.dTheta;
  rs = p.Results.rs;
  thetas = p.Results.thetas;
  maxRadius = p.Results.maxRadius;
  method = p.Results.method;
  extrapval = p.Results.extrapval;

  if numel( dr ) < 1, dr = 1; end;
  if numel( dTheta ) < 1, dTheta = defaultDTheta; end;
  if numel( maxRadius ) < 1, maxRadius=0; end;
  if numel( method ) < 1, method='linear'; end;
  if numel( extrapval ) < 1, extrapval = 0; end;

  sImg = size( img );
  if maxRadius==0, maxRadius = floor(min(size(img))/2); end;

  if numel( rs ) < 1, rs = 0:dr:maxRadius; end;
  if numel(thetas) < 1, thetas = -pi:dTheta:pi-dTheta; end;

  [rsImg,thetasImg] = ndgrid( rs, thetas );

  polarXs = rsImg .* cos(thetasImg);
  polarYs = rsImg .* sin(thetasImg);

  tmp = size2imgCoordinates( sImg );
  ys=tmp{1};  xs=tmp{2};
  [ys, xs] = ndgrid( ys, xs );

  outImg = interp2( xs, ys, img, polarXs, polarYs, method, extrapval );
end
