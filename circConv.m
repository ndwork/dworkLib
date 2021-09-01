
function out = circConv( inArray, kernel, op )
  % out = circConv( inArray, kernel [, op ] )
  %
  % Calculates the circular convolution of inArray with kernel.
  %   The origin of inArray is the first element.
  %   The origin of kernel is the center element.
  %
  % Inputs:
  % inArray - a 2D array
  % kernel - a small 2D array
  %
  % Optional Inputs:
  % op - either 'notransp' (default) or 'transp'
  %
  % Outputs:
  % out - a 2D array
  %
  % Written by Nicholas Dwork - Copyright 2021
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  if nargin < 3, op = 'notransp'; end

  sA = size( inArray );
  sK = size( kernel );
  newSize = max( sA, sK );

  pad1 = zeros( newSize );
  pad2 = zeros( newSize );
  str2eval1 = 'pad1( 1:sA(1)';
  str2eval2 = 'pad2( 1:sK(1)';
  for dim = 2 : numel( sA )
    str2eval1 = [ str2eval1, ', 1:sA(', num2str(dim), ')' ];   %#ok<AGROW>
    str2eval2 = [ str2eval2, ', 1:sK(', num2str(dim), ')' ];   %#ok<AGROW>
  end
  str2eval1 = [ str2eval1, ' ) = inArray;' ];
  str2eval2 = [ str2eval2, ' ) = kernel;' ];
  
  eval( str2eval1 );
  eval( str2eval2 );

  pad2 = circshift( pad2, -floor( sK / 2 ) );
  
  if strcmp( op, 'transp' )
    pad2 = flipAboutIndx( pad2, ones( 1, ndims( pad2 ) ) );
  end

  out = ifftn( fftn( pad1 ) .* fftn( pad2 ) );
end

