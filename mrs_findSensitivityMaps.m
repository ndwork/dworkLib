
function sMaps = mrs_findSensitivityMaps( coilSpectrums, noiseCoords )
  % sMaps = mrs_findSensitivityMaps( coilSpectrums )
  %
  % Inputs:
  % coilSpectrums - an array of size sImg x nCoils x nF representing the Fourier value
  %   of each spatial location and frequency bin
  %   The data is in space and temporal frequency
  %   sImg - a 1D array representing the size of the image
  %   nCoils - the number of coils
  %   nF - the number of frequency bins in the spectrum
  %
  % Written by Nicholas Dwork, Copyright 2021
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  sCoilSpectrums = size( coilSpectrums );
  sImg = sCoilSpectrums( 1 : end - 2 );
  nCoils = sCoilSpectrums( 3 );
  %nFreqs = sCoilSpectrums( 4 );

  noise = coilSpectrums( noiseCoords(1) : noiseCoords(2), noiseCoords(3) : noiseCoords(4), :, : );  
  noise = squeeze( sqrt( sum( abs( noise ).^2, 3 ) ) );

  nImg = prod( sImg );
  coilSpectrums = reshape( coilSpectrums, [ nImg sCoilSpectrums( end-1 : end ) ] );
  rsos = squeeze( sqrt( sum( abs( coilSpectrums ).^2, 2 ) ) );  % root-sum-of-squares
  denoised = removeRicianMean( rsos, noise );

  sMaps = zeros( [ nImg nCoils ] );
  for j = 1 : nImg
    %a = squeeze( sqrt( sum( coilSpectrums( j, :, : ).^2, 2 ) ) );
    a = transpose( denoised( j, : ) );
    b = transpose( squeeze( coilSpectrums( j, :, : ) ) );

    sMaps( j, : ) = a \ b;
  end

  sMaps = reshape( sMaps, [ sImg nCoils ] );
end

