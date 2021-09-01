
function out = crossCorrelate( inArray, kernel, op )
  % out = out = crossCorrelate( inArray, kernel [, op ] )
  %
  % Calculates the cross correlation of inArray with kernel.
  %   The origin of inArray is the first element.
  %   The origin of kernel is the center element.
  %
  % Inputs:
  % inArray - an N dimensional array
  % kernel - a small N dimensional array
  %
  % Optional Inputs:
  % op - either 'notransp' (default) or 'transp'
  %
  % Outputs:
  % out - a N dimensional array of size max( size( inArray ), size( kernel ) )
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
  newSize = max( sA + sK -1 );

  pad1 = zeros( newSize );
  pad2 = zeros( newSize );
  str2eval1 = 'pad1( 1:sA(1)';
  str2eval2 = 'pad2( 1:sK(1)';
  str2evalOut = 'out = out( 1:sA(1)';
  for dim = 2 : numel( sA )
    str2eval1 = [ str2eval1, ', 1:sA(', num2str(dim), ')' ];   %#ok<AGROW>
    str2eval2 = [ str2eval2, ', 1:sK(', num2str(dim), ')' ];   %#ok<AGROW>
    str2evalOut = [ str2evalOut, ', 1:sA(', num2str(dim), ')' ];  %#ok<AGROW>
  end
  str2eval1 = [ str2eval1, ' ) = inArray;' ];
  str2eval2 = [ str2eval2, ' ) = kernel;' ];
  str2evalOut = [ str2evalOut, ' );' ];

  eval( str2eval1 );
  eval( str2eval2 );

  pad2 = circshift( pad2, -floor( sK / 2 ) );
  
  if strcmp( op, 'notransp' )
    pad2 = flipAboutIndx( pad2, ones( ndims( inArray ), 1 ) );
  end

  out = ifftn( fftn( pad1 ) .* fftn( pad2 ) );
  eval( str2evalOut );
end

