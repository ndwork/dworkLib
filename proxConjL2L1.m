
function out = proxConjL2L1( x, t )
  % out = proxConjL2L1( x, t )
  %
  % Calculates the proximal operator of the conjugate of f(x) = t * L2L1( x )
  %
  % Inputs:
  % x - an N Dimensional array where the last dimension represents the groups
  % t - a scalar
  %
  % Written by Nicholas Dwork - Copyright 2020
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  nx = norms( x );
  scalings = ones( size( nx ) );
  tInv = 1 / t;
  scalings( nx > 1 ) = tInv ./ nx( nx > 1 );

  out = bsxfun( @times, x, scalings );
end
