
function [vShift, hShift, rotation] = findTransRotWithPCC( img1, img2, ...
  varargin )
  % [vShift, hShift, rotation] = findTransRotWithPCC( img1, img2 [, ...
  %   'dTheta', dTheta] )
  %
  % Finds translation and rotation so that
  %   img2 = Rotation( Translation( img1 ) )
  %
  % Output:
  % vShift = the vertical shift (in pixels)
  % hShift = the horizontal shift (in pixels)
  % rotation = the rotation (in radians)
  %
  % Written by Nicholas Dwork - Copyright 2017
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  defaultDTheta = [];
  p = inputParser;
  p.addParameter( 'dTheta', defaultDTheta, @isnumeric );
  p.parse( varargin{:} );
  dTheta = p.Results.dTheta;
  
  fftImg1 = fftshift(fft2( img1 ) );  absFft1 = abs( fftImg1 );
  fftImg2 = fftshift(fft2( img2 ) );  absFft2 = abs( fftImg2 );

  rotation = findRotWithPCC( absFft1, absFft2, 'dTheta', dTheta );
  rot2 = imrotate( img2, -rotation*180/pi, 'crop' );
  [vShift,hShift] = findTransWithPCC( img1, rot2 );
end