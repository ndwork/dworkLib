
function ihs = rgb2ihs( rgb )
  % ihs = rgb2ihs( rgb )
  % Converts an RGB color image to a Intensity-Hue-Saturation color image
  %
  % Written by Nicholas Dwork - Copyright 2018
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  sImg = size( rgb );

  nPixels = prod( sImg(1:2) );
  rgbVec = zeros( 3, nPixels );
  rgbVec(1,:) = reshape( rgb(:,:,1), [1 nPixels] );
  rgbVec(2,:) = reshape( rgb(:,:,2), [1 nPixels] );
  rgbVec(3,:) = reshape( rgb(:,:,3), [1 nPixels] );

  M = [ 1/sqrt(3) 1/sqrt(3) 1/sqrt(3); ...
        1/sqrt(6) 1/sqrt(6) -2/sqrt(6); ...
        1/sqrt(2) -1/sqrt(2) 0; ];
  ihsVec = M * rgbVec;

  ihs = zeros( size(rgb) );
  ihs(:,:,1) = reshape( ihsVec(1,:), sImg(1:2) );
  ihs(:,:,2) = reshape( ihsVec(2,:), sImg(1:2) );
  ihs(:,:,3) = reshape( ihsVec(3,:), sImg(1:2) );

end

