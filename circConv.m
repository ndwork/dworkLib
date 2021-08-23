
function out = circConv( inArray, kernel )
  % out = circConv( inArray, kernel )
  %
  % Calculates the circular convolution of inArray with kernel.
  %   The origin of inArray is the first element.
  %   The origin of kernel is the center element.
  %
  % Inputs:
  % inArray - a 2D array
  % kernel - a small 2D array
  %
  % Outputs:
  % out - a 2D array that is the circular convolution of inArray with the kernel
  %
  % Written by Nicholas Dwork - Copyright 2021
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  s1 = size( inArray );
  s2 = size( kernel );
  newSize = max( s1, s2 );

  pad1 = zeros( newSize );
  pad2 = zeros( newSize );
  str2eval1 = 'pad1( 1:s1(1)';
  str2eval2 = 'pad2( 1:s2(1)';
  for dim = 2 : numel( s1 )
    str2eval1 = [ str2eval1, ', 1:s1(', num2str(dim), ')' ];   %#ok<AGROW>
    str2eval2 = [ str2eval2, ', 1:s2(', num2str(dim), ')' ];   %#ok<AGROW>
  end
  str2eval1 = [ str2eval1, ' ) = inArray;' ];
  str2eval2 = [ str2eval2, ' ) = kernel;' ];
  
  eval( str2eval1 );
  eval( str2eval2 );

  pad2 = circshift( pad2, -floor( s2 / 2 ) );

  out = ifftn( fftn( pad1 ) .* fftn( pad2 ) );
end

