
function [FVals,phantImg] = makeNickPhantom( traj )
  % [FVals,phantImg] = makeNickPhantom( traj )
  %
  % Inputs:
  % traj - a Kx2 array, where K is the number of trajectory points
  %   the first/second column is ky / kx
  %
  % Outputs:
  % FVals - An array of size Kx1 with Fourier values
  % phantImg - an image of size [256 256] that is the image
  %
  % Written by Nicholas - Copyright 2017
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp( 'Usage:  [FVals,phantImg] = makeNickPhantom( traj )' );
    if nargout > 0, FVals = []; end
    if nargout > 1, phantImg = []; end
    return
  end

  %-- Phantom parameters
  M = 256;
  N = 256;

  xRect = 35;     widthRect = 120;
  yRect = 0;      heightRect = 75;
  valRect = 1;

  xCirc = -45;    circRad = 50;
  yCirc = 45;     circVal = 2.0;

  xTri = -50;     widthTri = 70;
  yTri = -50;     heightTri = 70;
  triMaxVal = 3.0;

  xRect_2 = 70;  widthRect_2 = 20;
  yRect_2 = 0;    heightRect_2 = 200;
  valRect_2 = 2;

  %-- Compute the Fourier transform values
  FValsRect = valRect * exp( -1i * 2*pi * ( xRect * traj(:,2) + yRect * traj(:,1) ) ) .* ...
    (widthRect * heightRect) .* sinc( widthRect * traj(:,2) ) .* sinc( heightRect * traj(:,1) );

  circDiam = 2 * circRad;
  FValsCirc = circVal * exp( -1i * 2*pi * ( xCirc * traj(:,2) + yCirc * traj(:,1) ) ) .* ...
    (circDiam * circDiam) .* jinc( circDiam * sqrt(traj(:,2).^2 + traj(:,1).^2) );

  FValsTri = triMaxVal * exp( -1i * 2*pi * ( xTri * traj(:,2) + yTri * traj(:,1) ) ) .* ...
    (widthTri * heightTri) .* sinc(widthTri * traj(:,2)).^2 .* sinc( heightTri * traj(:,1) ).^2;

  FValsRect_2 = valRect_2 * exp( -1i * 2*pi * ( xRect_2 * traj(:,2) + yRect_2 * traj(:,1) ) ) .* ...
    (widthRect_2 * heightRect_2) .* sinc( widthRect_2 * traj(:,2) ) .* sinc( heightRect_2 * traj(:,1) );

  FVals = FValsRect + FValsCirc + FValsTri + FValsRect_2;


  %-- Create an image of the phantom
  if nargout > 1
    imgXs = ones(M,1) * (1:N) - (floor(N/2)+1);
    imgYs = (1:M)' * ones(1,N) - (floor(M/2)+1);

    xImgRect = abs( (imgXs - xRect) / widthRect ) < 0.5;
    yImgRect = abs( (imgYs - yRect) / heightRect ) < 0.5;
    imgRect = valRect * xImgRect .* yImgRect;

    imgCirc = sqrt( (imgXs-xCirc).^2 + (imgYs-yCirc).^2 ) < circRad;
    imgCirc = imgCirc * circVal;

    imgTriX = max( 1 - abs( (imgXs-xTri) / widthTri), 0 );
    imgTriY = max( 1 - abs( (imgYs-yTri) / heightTri), 0 );
    imgTri = triMaxVal * imgTriX .* imgTriY;

    xImgRect_2 = abs( ( imgXs - xRect_2 ) / widthRect_2 ) < 0.5;
    yImgRect_2 = abs( ( imgYs - yRect_2 ) / heightRect_2 ) < 0.5;
    imgRect_2 = valRect_2 * xImgRect_2 .* yImgRect_2;

    phantImg = imgRect + imgCirc + imgTri + imgRect_2;

  end

end

