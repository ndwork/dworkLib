
function [ fsr, sFSR ] = mri_makeFullySampledCenterRegion( sImg, wavSplit )
  % [ fsr, sFSR ] = mri_makeFullySampledCenterRegion( sImg, wavSplit )
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It 
  % is offered without any warranty expressed or implied, including the 
  % implied warranties of merchantability or fitness for a particular 
  % purpose.

  wavMask = makeLowFreqWavMask( sImg, wavSplit );

  lastX = find( wavMask(1,:), 1, 'last' );
  lastY = find( wavMask(:,1), 1, 'last' );

  fsr = circshift( wavMask, [ -floor(lastY * 0.5) -floor(lastX * 0.5) ] );
  fsr = fftshift( fsr );

  if nargout > 1
    sFSR = sImg ./ ( 2.^( log( size(wavSplit) * 2 ) ./ log(2) ) );
  end

end

