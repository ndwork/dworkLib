
function out = proxL2Sq( v, t, b, A )
  % out = proxL2Sq( v, t, b, A )
  %
  % Evaluate the proximal operator of tf where f(x) = 1/2 || A x - b ||_2^2
  %
  % Inputs:
  % x - the argument of the proximal operator
  %     Note: could be a scalar or multi-dimensional array
  %
  % Optional Inputs:
  % t - scaling of the L2 norm squared function
  %     Either a scalar or a vector that is the same size as v
  % b - the shifting factor
  % A - the matrix
  %
  % Output:
  % out - the result of the proximal operator
  %
  % Written by Nicholas Dwork, Copyright 2019
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  switch nargin
    case 1
      out = 0.5 .* v;
    case 2
      out = v ./ ( 1 + t );
    case 3
      out = ( v + t .* b ) ./ ( 1 + t );
    case 4
      if isa( A, 'function_handle' )
        % solve with the conjugate gradient method
        E = @(x) A( A( x ), 'transp' ) + ( ones( size(x) ) ./ t );
        v = A( b, 'transp' ) + ( v ./ t );
        out = cgs( E, v );
        % TODO: If A is tight frame, this can be solved much faster
      else
        out = ( A'*A + eye( size(A,2) ) ./ t ) \ ( A' * b + v ./ t );
      end
  end

end
