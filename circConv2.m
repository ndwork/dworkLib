
function out = circConv2( A, K, op )
  % out = circConv2( A, K [, op ] )
  %
  % Calculates the two-dimensional circular convolution of A and K
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
    disp( 'Usage:  out = circConv2( A, K [, op, ''notransp''/''transp'' ] )' );
    if nargout > 0, out = []; end
    return
  end

  sA = size( A );
  sK = size( K );

  sOut = max( sA, sK );
  padA = padData( A, sOut );
  padK = padData( K, sOut );

  if nargin > 2 && strcmp( op, 'transp' )
    padA = flip( flip( conj( padA ), 1 ), 2 );
    if mod( size( padA, 1 ), 2 ) == 0
      padA = circshift( padA, 1, 1 );
    end
    if mod( size( padA, 2 ), 2 ) == 0
      padA = circshift( padA, 1, 2 );
    end
  end

  out = fftshift2( ifft2( fft2( ifftshift2( padA ) ) .* fft2( ifftshift2( padK ) ) ) );
end

