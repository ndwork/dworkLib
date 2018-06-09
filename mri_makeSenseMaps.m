
function maps = mri_makeSenseMaps( kData )
  % Compute sensitivity maps as described in "SENSE: Sensitivity Encoding
  % for Fast MRI" by Pruessmann et al. (1999)
  %
  % Inputs:
  % kData is a 4D complex array of size (Ny, Nx, nSlices, nCoils )
  %
  % Outputs:
  % maps is a 4D complex array: (Ny, Nx, nSlices, nCoils )
  %
  % Written by Nicholas Dwork - Copyright 2018
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  [Ny,Nx,nSlices,nCoils] = size( kData );

  % Reconstruct the images
  coilRecons = zeros( Ny, Nx, nSlices, nCoils );
  ssqRecon = zeros( Ny, Nx, nSlices );
  for slice=1:nSlices

    tmpReconSq = zeros( Ny, Nx );
    for coil=1:nCoils
      tmpFFT = squeeze( kData(:,:,slice,coil) );
      tmpImg = fftshift( ifft2( ifftshift(tmpFFT) ) );

      coilRecons(:,:,slice,coil) = tmpImg;
      tmpReconSq = tmpReconSq + conj(tmpImg) .* tmpImg;
    end

    ssqRecon(:,:,slice) = sqrt( tmpReconSq );
  end

  % Determine the senstivity maps
  maps = zeros( Ny, Nx, nSlices, nCoils );
  for coil = 1:nCoils
    maps(:,:,:,:,coil) = coilRecons(:,:,:,:,coil) ./ ssqRecon;
  end

end
