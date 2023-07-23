
function recon = mri_reconSSQ( kData, varargin )
  % recon = mri_reconSSQ( kData [, 'multiSlice', true/false ] )
  %
  % Perform a sum of squared reconstruction
  %
  % Inputs:
  % kData is an array of size ( Ny, Nx, nSlices, ..., nCoils ) of kSpace values
  %
  % Optional Inputs:
  % multiSlice - if set to true, assumed that each slice contains its own Fourier data
  %   Otherwise, assumed it's a 3D reconstruction (and a Fourier transform accross 3rd dimension
  %   is necessary).
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

  if nargin < 1
    disp( 'Usage:  recon = mri_reconSSQ( kData [, ''multiSlice'', true/false ] )' );
    return;
  end

  coilRecons = mri_reconIFFT( kData, varargin{:} );

  tmpReconSq = coilRecons .* conj( coilRecons );

  recon = sqrt( sum( tmpReconSq, ndims(kData) ) );

end
