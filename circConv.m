
function out = circConv( A, K, op )
  % out = circConv( A, K [, op, 'notransp/transp' ] )
  %
  % Calculates the circular convolution of inArray with kernel.
  %   The origin of inArray is the first element.
  %   The origin of kernel is the center element.
  %
  % Optional Inputs:
  % op - either 'notransp' (default) or 'transp'
  %         If 'transp' is set, calculates the adjoint of f(x) = circConv( A, x )
  %
  % Inputs:
  % A - a (possibly) multidimensional array
  % K - a (possibly) multidimensional array of the same number of dimensions as a
  %     sometimes referred to as the convolution kernel
  %
  % Outputs:
  % out - an array of size equal to the maximum size of a and b in each dimension that
  %       represents the result of a circular convolution  
  %
  % Written by Nicholas Dwork - Copyright 2021
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  if nargin < 1
    disp( 'Usage:  out = circConv( A, K [, op, ''notransp/transp'' )' );
    if nargout > 0, out = []; end
    return
  end

  sA = size( A );
  sK = size( K );

  sOut = max( sA, sK );
  padA = padData( A, sOut );
  padK = padData( K, sOut );

  if nargin > 2  &&  strcmp( op, 'transp' )
    padA = flipAllDims( conj( padA ) );
    padA = circshift( padA, int8( mod( sOut, 2 ) == 0 ) );
  end

  out = fftshift( ifftn( fftn( ifftshift( padA ) ) .* fftn( ifftshift( padK ) ) ) );
end

