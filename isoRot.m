
function out = isoRot( img, theta )
  % out = isoRot( img, theta )
  %
  % Performs a rotation that is an isometry by approximating a rotation as a set of shears.
  % as described here:
  %   http://www.numerical-tours.com/matlab/inverse_8_tomography/
  % Performing isoRot with negative theta is the invers of isoRot with theta.
  % This rotation is a true isometry and is invertible.
  %
  % Written by Nicholas Dwork - Copyright 2013
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp( 'Usage:  out = isoRot( img, theta )' );
    return
  end

  absTheta = abs( theta );
  nPiHalves = floor( absTheta / ( 0.5 * pi ) );
  newTheta = absTheta - ( nPiHalves * 0.5 * pi );
  nPiEighths = floor( newTheta / ( 0.125 * pi ) );
  extraTheta = newTheta - nPiEighths * 0.125 * pi;

  if theta >= 0    
    out = rot90( img, nPiHalves );
    for i = nPiEighths
      out = isoRot_theta( out, 0.125 * pi );
    end
    out = isoRot_theta( out, extraTheta );

  else  % theta < 0
    out = isoRot_theta( img, -extraTheta );
    for i = nPiEighths
      out = isoRot_theta( out, -0.125 * pi );
    end
    out = rot90( out, -nPiHalves );

  end
end

function out = isoRot_theta( img, theta )
  sh = shearX( img, -tan(theta/2) );
  sh = shearY( sh, sin(theta) );
  out = shearX( sh, -tan(theta/2) );
end

function sheared = shearX( img, a )
  % sheared = shearX( f, a )
  % A shear as implemented here:
  % https://www.ceremade.dauphine.fr/~peyre/numerical-tour/tours/inverse_8_tomography/

  sh = shearY( img', a );
  sheared = sh';
end

function sheared = shearY( img, a )
  [Ny, Nx] = size(img);

  tx = [0:Nx/2-1 0 -Nx/2+1:-1]';
  ty = [0:Ny/2-1 0 -Ny/2+1:-1]';
  [X,Y] = meshgrid(tx,ty);

  shiftedF = fftshift(img);
  tmp = real( ifft( fft(shiftedF) .* exp(a*2i*pi*X.*Y/Ny) ) );
  sheared = ifftshift( tmp );
end
