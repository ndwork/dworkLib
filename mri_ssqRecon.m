
function recon = mri_ssqRecon( kData )
  % Perform a sum of squared reconstruction
  %
  % Inputs:
  % kData is an array of size ( Ny, Nx, nSlices, ..., nCoils ) of kSpace values
  %
  % Output:
  % recon is the reconstructed image
  %
  % Written by Nicholas Dwork - Copyright 2018
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  coilRecons = mri_fftRecon( kData );

  tmpReconSq = coilRecons .* conj( coilRecons );
  
  recon = sqrt( sum( tmpReconSq, ndims(kData) ) );

end
