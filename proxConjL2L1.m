
function out = proxConjL2L1( x, ~ )
  % out = proxConjL2L1( x, ~ )
  %
  % Calculates the proximal operator of the L2,L1 norm applied to the input
  %
  % Written by Nicholas Dwork - Copyright 2020
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  nx = norms( x, 2 );
  scalings = ones( size( nx ) );
  scalings( nx > 1 ) = 1 ./ scalings( nx > 1 );

  out = bsxfun( @times, x, scalings );
end
