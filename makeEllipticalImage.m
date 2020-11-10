
function out = makeEllipticalImage( sImg, varargin )
  % out = makeEllipticalImage( sImg, [ rx, ry, rot ] )
  %
  % Creates an elliptical image.  The zero-levelset is where 
  % x^2/rx + y^2/ry = 1
  %
  % Inputs:
  % sImg - two element array specifying the size of the image [Ny,Nx]
  % rx - the length of the semimajor axis in the horizontal direction
  % ry - the length of the semimajor axis in the vertical direction
  % rot (optional) - clockwise rotation in radians
  %
  % Written by Nicholas Dwork
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.


  defaultRot = 0;
  defaultRx = floor(sImg(2)/2);
  defaultRy = floor(sImg(1)/2);
  p = inputParser;
  p.addOptional( 'rx', defaultRx, @isnumeric );
  p.addOptional( 'ry', defaultRy, @isnumeric );
  p.addOptional( 'rot', defaultRot, @isnumeric );
  p.parse( varargin{:} );
  rx = p.Results.rx;
  ry = p.Results.ry;
  rot = p.Results.rot;

  y = (1:sImg(1)) - ceil((sImg(1)+1)/2);
  x = (1:sImg(2)) - ceil((sImg(2)+1)/2);
  [xs,ys] = meshgrid( x, y );

  if rot~=0
    pts = [ transpose(xs(:)); transpose(ys(:)); ];
    R = [ cos(rot) sin(rot); -sin(rot), cos(rot); ];
    Rpts = R * pts;
    xs = Rpts(1,:);  xs=reshape( xs, sImg );
    ys = Rpts(2,:);  ys=reshape( ys, sImg );
  end

  out = ( (xs.*xs)/(rx*rx) + (ys.*ys)/(ry*ry) ) - 1;
end
