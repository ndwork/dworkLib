
function v = fitPolyToData( N, x, y )
  % v = fitPolyToData( N, x, y ) or
  % v = fitPolyToData( N, y )
  %
  % This function finds a polynomial p so that ||y-p(x)||_2 is minimized
  %
  % Inputs:
  % N - the order of the polynomial
  % x - (optional) domain values
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

  if nargin < 2
    y = x;
    x = 1:numel(y);
  end

  A = ones( numel(x), N+1 );
  for i=1:N
    A(:,i+1) = A(:,i) .* x(:);
  end

  v = A \ y(:);
end
