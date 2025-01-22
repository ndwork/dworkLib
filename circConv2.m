
function out = circConv2( A, K, varargin )
  % out = circConv2( A, K [, op, 'ndimsOut', ndimsOut ] )
  %
  % Calculates the two-dimensional circular convolution of A and K
  %
  % Optional Inputs:
  % op - either 'notransp' (default) or 'transp'
  % ndimsOut - the number of dimensions of the output, used when 'transp' is set
  %   This is used when the inputs of the forward operator have different dimensions
  %   By default, assumes the same dimensions as max( ndims(A), ndims(K) )
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

  if nargin < 2
    disp( 'Usage:  out = circConv2( A, K [, op, ''notransp''/''transp'' ] )' );
    if nargout > 0, out = []; end
    return
  end

  p = inputParser;
  p.addOptional( 'op', 'notransp', @(x) true );
  p.addParameter( 'ndimsOut', [], @ispositive );
  p.parse( varargin{:} );
  op = p.Results.op;
  ndimsOut = p.Results.ndimsOut;

  sA = size( A );
  sK = size( K );

  sOut = max( sA(1:2), sK(1:2) );
  if any( sA(1:2) ~= sOut )
    sA_new = [ sOut sA(3:end) ];
    A = padData( A, sA_new );
  end
  if any( sK(1:2) ~= sOut )
    sK_new = [ sOut sK(3:end) ];
    K = padData( K, sK_new );
  end

  if nargin > 2 && strcmp( op, 'transp' )
    A = flipDims( conj( A ), 'dims', [ 1 2 ] );
  end

  %out = fftshift2( ifft2( fft2( ifftshift2( A ) ) .* fft2( ifftshift2( K ) ) ) );
  out = fftshift2( ifft2( bsxfun( @times, fft2( ifftshift2( A ) ), fft2( ifftshift2( K ) ) ) ) );

  if nargin > 2 && strcmp( op, 'transp' ) && numel( ndimsOut ) > 0
    out = sum( out, ndimsOut+1 : ndims(out) );
  end
end

