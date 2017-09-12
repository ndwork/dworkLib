
function rotation = findRotWithPCC( img1, img2, varargin )
  % rotation = findRotWithPCC( img1, img2 [, 'dTheta', dTheta] )
  % Note: rotation must be -pi/2 and pi/2
  %
  % Inputs:
  % dTheta - resolution of rotation angle  (default is 1/180*pi radians)
  %
  % Output:
  % rotation - counterclockwise rotation to apply to img1 (in radians) so
  %            that it is aligned with img2
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  defaultDTheta = 1*pi/180;
  p = inputParser;
  p.addParameter( 'dTheta', defaultDTheta, @isnumeric );
  p.parse( varargin{:} );
  dTheta = p.Results.dTheta;
  
  if numel( dTheta ) == 0, dTheta = defaultDTheta; end;
  
  lp1 = smoothImg( img1, 49, 'gaussian', 19 );  hp1 = img1-lp1;
  lp2 = smoothImg( img2, 49, 'gaussian', 19 );  hp2 = img2-lp2;

  thetas = 0:dTheta:pi-dTheta;
  fftHp1 = fftshift( fft2( hp1 ) );
  fftHp2 = fftshift( fft2( hp2 ) );
  
  polar1 = img2Polar( abs( fftHp1 ), 'thetas', thetas );
  polar2 = img2Polar( abs( fftHp2 ), 'thetas', thetas );

  pcc = phaseCrossCorrelate( polar2, polar1 );
  [~,maxIndx] = max( pcc(1,:) );
  
  rotation = thetas( maxIndx );
  if rotation > pi/2, rotation = rotation-pi; end;
end
