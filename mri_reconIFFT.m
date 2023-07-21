
function recons = mri_reconIFFT( kData, varargin )
  % recons = mri_reconIFFT( kData [, 'multiSlice', true/false ] )
  %
  % Performs an inverse FFT of each coil
  %
  % Inputs:
  % kData is an array of size ( Ny, Nx, nSlices, ..., nCoils )
  %
  % Output:
  % recons is the reconstructed image, an array of the same size as the input
  %
  % Optional Inputs:
  % multiSlice - if set to true, assumed that each slice contains its own Fourier data
  %   Otherwise, assumed it's a 3D reconstruction (and a Fourier transform accross 3rd dimension
  %   is necessary).
  %
  % Written by Nicholas Dwork - Copyright 2018
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp('Usage:  recons = mri_reconIFFT( kData [, ''multiSlice'', true/false );' );
    if nargout > 0, recons = []; end
    return;
  end

  p = inputParser;
  p.addParameter( 'multiSlice', false, @islogical );
  p.parse( varargin{:} );
  multiSlice = p.Results.multiSlice;

  recons = fftshift2( ifft2( ifftshift2( kData ) ) );

  if ~ismatrix( kData ) && multiSlice == false
    recons = fftshift( ifft( ifftshift( recons, 3 ), [], 3 ), 3 );
  end

end
