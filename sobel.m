
function out = sobel( img )
  % out = sobel( img )
  %
  % Calculates the sobel operator of an image
  %
  % Written by Nicholas Dwork - Copyright 2020
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp( 'Usage: out = sobel( img )' );
  end

  Gx = [ [ 1 0 -1 ]; ...
         [ 2 0 -2 ]; ...
         [ 1 0 -1 ]; ];
  sxImg = conv2( img, Gx, 'same' );

  Gy = [ [  1  2  1 ]; ...
         [  0  0  0 ]; ...
         [ -1 -2 -1 ]; ];
	syImg = conv2( img, Gy, 'same' );

  out = sqrt( sxImg.^2 + syImg.^2 );
end
