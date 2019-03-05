
function fused = hsiFusion( colorImg, monoImg )
  % fused = ihsFusion( colorImg, monoImg )
  %
  % performs the fast IHS fusion algorithm of "A New Intensity-Hue-Saturation 
  %   Fusion Approach to Image Fusion With a Tradeoff Parameter" by Choi,
  %   2006
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  Intensity = mean( colorImg, 3 );

  fused = zeros( size(colorImg) );

  diffMI = monoImg - Intensity;
  fused(:,:,1) = colorImg(:,:,1) + diffMI;
  fused(:,:,2) = colorImg(:,:,2) + diffMI;
  fused(:,:,3) = colorImg(:,:,3) + diffMI;

  fused = min( max( fused, 0 ), 1 );
end
