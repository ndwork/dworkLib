
function out = makeEllipticalImage( sImg, ry, rx, varargin )
  % out = makeEllipticalImage( sImg, ry, rx [, rot ] )
  %
  % Creates an elliptical image.  The zero-levelset is where 
  % x^2/rx + y^2/ry = 1
  %
  % Inputs:
  % sImg - two element array specifying the size of the image [Ny,Nx]
  % ry - the length of the semimajor axis in the y direction
  % rz - the length of the semimajor axis in the x direction
  % rot (optional) - clockwise rotation in radians
  %
  % Written by Nicholas Dwork
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.


  defaultRot = 0;
  p = inputParser;
  p.addOptional( 'rot', defaultRot );
  p.parse( varargin{:} );
  rot = p.Results.rot;

  y = (1:sImg(1)) - ceil((sImg(1)+1)/2);
  x = (1:sImg(2)) - ceil((sImg(2)+1)/2);
  [xs,ys] = meshgrid( x, y );
  
  if rot~=0,
    pts = [ transpose(xs(:)); transpose(ys(:)); ];
    R = [ cos(rot) sin(rot); -sin(rot), cos(rot); ];
    Rpts = R * pts;
    xs = Rpts(1,:);  xs=reshape( xs, sImg );
    ys = Rpts(2,:);  ys=reshape( ys, sImg );
  end

  out = ( (xs.*xs)/(rx*rx) + (ys.*ys)/(ry*ry) ) - 1;
end
