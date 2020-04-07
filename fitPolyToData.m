
function [v,RSquared] = fitPolyToData( N, x, y )
  % [v,RSquared] = fitPolyToData( N, x [, y] )
  %
  % This function finds a polynomial p so that || y - p(x) ||_2 is minimized
  %
  % Inputs:
  % N - the order of the polynomial
  % x - (optional) domain values.  If not supplied, x = 1, 2, ..., numel(y)
  % y - range values
  %
  % Outputs:
  % v = a vector of polynomial coefficients.
  %   p(x) = v(1) + v(2)*x + v(3)*x^3 + ... + v(N)*x^(N)
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp( 'Usage:  [v,RSquared] = fitPolyToData( N, x [, y] )' );  return;
  end
  if nargout > 1 && N ~= 1
    error( 'N must be 1 to return RSquared' );  return;
  end

  if nargin < 3
    y = x;
    x = 1:numel(y);
  end

  A = ones( numel(x), N+1 );
  for i=1:N
    A(:,i+1) = A(:,i) .* x(:);
  end

  v = A \ y(:);
  
  if nargout > 1
    p = evaluatePoly( v, x );
    meanY = mean( y );
    ssTot = sum( ( y(:) - meanY ).^2 );
    ssReg = sum( ( p(:) - meanY ).^2 );
    RSquared = ssReg / ssTot;
  end
end
