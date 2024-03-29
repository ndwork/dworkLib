
function out = proxCompositionAffine( proxf, x, A, b, alpha, t )
  % Computes the proximal operator of t * f( A x + b ) where A A^T = (1/alpha) I.
  %
  % Inputs:
  % proxf - a function handle to the proximal operator of f
  % x - the domain value for proximal operator evaluation
  % A - either a matrix or a function handle that applies the matris
  % b - the shifting of the domain variable
  % alpha - a scalar
  % t - a scalar
  %
  % Written by Nicholas Dwork, Copyright 2024
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if isa( A, 'function_handle' )
    Ax = A( x );
    out = x - alpha * A( Ax + b - proxf( Ax + b, t/alpha ) , 'transp' );
  else
    Ax = A * x;
    out = x - alpha * A' * ( Ax + b - proxf( Ax + b, t/alpha ) );
  end
end
