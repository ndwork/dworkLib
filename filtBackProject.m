
function out = filtBackProject( sino, sImg, varargin )
  % out = filtBackProject( sino, sImg [, 'detCenter', detCenter, 'dProjAngle', dProjAngle, ...
  %                        'dSize', dSize, 'projAngles', projAngles ] )
  %
  % Inputs:
  % sino - a parallel beam sinogram with one column per detector and one row per projection
  % sImg - the size of the output image
  %
  % Optional Inputs:
  % detCenter - the location of the center detector (defaults to center of projection)
  %   Note, for optical projection tomography, this is a component of the principal point
  % dProjAngle - the angle between adjacent rotation angles
  %   Either dAngle or projAngles must be supplied
  % dSize - the size of each detector.  (More specifically, the spacing between detectors.)
  % projAngles - the angle of the principal ray for each projection
  %   Either dAngle or projAngles must be supplied
  %
  % Written by Nicholas Dwork, Copyright 2025
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  % TODO: Density compensation may be beneficial

  p = inputParser;
  p.addParameter( 'detCenter', [] );
  p.addParameter( 'dProjAngle', [] );
  p.addParameter( 'dSize', 1, @ispositive );
  p.addParameter( 'projAngles', [] );
  p.parse( varargin{:} );
  detCenter = p.Results.detCenter;
  dProjAngle = p.Results.dProjAngle;
  dSize = p.Results.dSize;
  projAngles = p.Results.projAngles;

  if numel( dProjAngle ) == 0  &&  numel( projAngles ) == 0
    error( 'Either dProjAngle or projAngles must be provided.' );
  end

  nDetectors = size( sino, 1 );
  nProjections = size( sino, 2 );

  if numel( projAngles ) == 0
    projAngles = dProjAngle * ( 0 : nProjections-1 );
  end

  if numel( detCenter ) == 0
    n = size2imgCoordinates( nDetectors );
  else
    n = ( (0:nDetectors-1) - detCenter );
  end

  h = zeros( numel(n), 1 );
  oddNs = find( mod(n,2)~=0 );
  h(oddNs) = -1 ./ ( n(oddNs).*n(oddNs) * pi*pi * dSize*dSize );
  h( n==0 ) = 1 / ( 4 * dSize*dSize );
  nPadded = 2 * nDetectors;
  hZP = padData( h, nPadded );

  fftHZP = fft( ifftshift( hZP ) );
  fftHZP = fftHZP * ones( 1, nProjections );

  sinoZP = padData( sino, [ nPadded, nProjections ] );
  fftSino = fft( ifftshift( sinoZP, 1 ), [], 1 );

  fftFiltSino = fftHZP .* fftSino;
  filtSino = dSize * fftshift( ifft( fftFiltSino, [], 1 ), 1 );
  filtSino = cropData( filtSino, [ nDetectors, nProjections ] );

  out = backProject( filtSino, projAngles, sImg );
  out = out * pi / nProjections;
end
