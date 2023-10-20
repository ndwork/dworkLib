
function recons = mri_reconIFFT( kData, varargin )
  % recons = mri_reconIFFT( kData [, 'multiSlice', true/false, 'dims', dims ] )
  %
  % Performs an inverse FFT of each coil
  %
  % Inputs:
  % kData is an array of size ( Ny, Nx, nSlices, ..., nCoils )
  %
  % Optional Inputs:
  % dims - an array of the indices over which the IFFT should be taken
  % multiSlice - if set to true, then only does an IFFT on the first two dimensions
  %
  % Output:
  % recons is the reconstructed image, an array of the same size as the input
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
  p.addParameter( 'dims', [], @isnumeric );
  p.addParameter( 'multiSlice', false, @islogical );
  p.parse( varargin{:} );
  dims = p.Results.dims;
  multiSlice = p.Results.multiSlice;

  if numel( dims ) > 0
    recons = kData;
    for dimIndx = 1 : numel( dims )
      dim = dims( dimIndx );
      recons = fftshift( ifft( ifftshift( recons, dim ), [], dim ), dim );
    end
  else

    recons = fftshift2( ifft2( ifftshift2( kData ) ) );
  
    if ~ismatrix( kData ) && multiSlice == false
      recons = fftshift( ifft( ifftshift( recons, 3 ), [], 3 ), 3 );
    end
  end

end
