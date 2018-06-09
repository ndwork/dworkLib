
function imgData = mri_ssqRecon( kData )
  % Perform a sum of squared reconstruction
  %
  % Inputs:
  % kData is an array of size (Ny,Nx,nSlices,nCoils) of kSpace values
  %
  % Output:
  % recon is the reconstructed image
  %
  % Written by Nicholas Dwork - Copyright 2018
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  [Ny,Nx,nSlices,nCoils] = size( kData );
  imgData = zeros( Ny, Nx, nSlices );
  for slice=1:nSlices
    tmpReconSq = zeros( Ny, Nx );
    for coil=1:nCoils
      tmpFFT = squeeze( kData(:,:,slice,phase,coil) );
      tmpImg = fftshift( ifft2( ifftshift(tmpFFT) ) );
      tmpReconSq = tmpReconSq + conj(tmpImg) .* tmpImg;
    end
    imgData(:,:,slice,phase) = sqrt( tmpReconSq );
  end
end
