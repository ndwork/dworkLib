
function out = proxCompositionAffine( proxf, x, A, b, alpha, t )
  % Computes the proximal operator of t * f( A x + b ) where A A^T = (1/alpha) I.
  %
  % Inputs:
  % proxf - a function handle to the proximal operator of f
  % x - the domain value for proximal operator evaluation
  % A - either a matrix or a function handle that applies the matrix
  %     (default is identity)
  % b - the shifting of the domain variable
  %     (default is 0)
  % alpha - a scalar
  % t - a scalar
  %
  % Written by Nicholas Dwork, Copyright 2024
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if numel( b ) == 0, b = 0; end
  if numel( alpha ) == 0, alpha = 1; end
  if numel( t ) == 0, t = 1; end

  if numel( A ) > 0

    if isa( A, 'function_handle' )
      Ax = A( x );
      out = x - alpha * A( Ax + b - proxf( Ax + b, t/alpha ) , 'transp' );
    else
      Ax = A * x;
      out = x - alpha * A' * ( Ax + b - proxf( Ax + b, t/alpha ) );
    end

  else
    if alpha ~= 1, error( 'alpha must be 1 when A is the identity' ); end
    out = proxf( x + b, t ) - b;
  end
end
