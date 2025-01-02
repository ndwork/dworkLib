
function out  = circConv1( A, K, op )
  % out = circConv1( A, K [, op ] )
  %
  % Calculates the one-dimensional circular convolution of A and K
  %
  % Optional Inputs:
  % op - either 'notransp' (default) or 'transp'
  %
  % Inputs:
  % A - a (possibly) multidimensional array
  % K - a (possibly) multidimensional array of the same number of dimensions as a
  %     sometimes referred to as the convolution kernel
  %
  % Outputs:
  % out - an array of size equal to the maximum size of A and K in each dimension that
  %       represents the result of a circular convolution  
  %
  % Written by Nicholas Dwork - Copyright 2024
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  if nargin < 1
    disp( 'Usage:  out = circConv1( A, K [, op, ''notransp''/''transp'' ] )' );
    if nargout > 0, out = []; end
    return
  end

  nA = numel( A );
  nK = numel( K );

  nOut = max( nA, nK );
  padA = padData( A(:), nOut );
  padK = padData( K(:), nOut );

  if nargin > 2 && strcmp( op, 'transp' )
    if mod( nOut, 2 ) == 0
      padA = circshift( flip( conj( padA ) ), 1 );
    else
      padA = flip( conj( padA ) );
    end
  end

  out = fftshift( ifft( fft( ifftshift( padA ) ) .* fft( ifftshift( padK ) ) ) );
end
