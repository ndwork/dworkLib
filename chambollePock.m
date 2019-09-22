
function [xStar,objValues] = chambollePock( x, proxf, proxgConj, sigma, tau, varargin )
  % [xStar,objValues] = chambollePock( x, proxf, proxgConj, sigma, tau [, ...
  %   'A', A, 'f', f, 'g', g ] )
  %
  % minimizes f( x ) + g( A x )
  %
  % Optional Inputs:
  % A - if A is not provided, it is assumed to be the identity
  % f - to determine the objective values, f must be provided
  % g - to determine the objective values, g must be provided
  % N - the number of iterations that ADMM will perform (default is 100)
  %
  % Outputs:
  % xStar - the optimal point
  %
  % Optional Outputs:
  % objValues - a 1D array containing the objective value of each iteration
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addParameter( 'A', [] );
  p.addParameter( 'f', [] );
  p.addParameter( 'g', [] );
  p.addParameter( 'N', 100, @ispositive );
  p.addParameter( 'theta', 1, @(x) x > 0 && x < 2 );
  p.parse( varargin{:} );
  A = p.Results.A;
  f = p.Results.f;
  g = p.Results.g;
  N = p.Results.N;

  if numel( A ) == 0
    applyA = @(x) x;
    applyAT = @(x) x;
  elseif isnumeric( A )
    applyA = @(x) A * x;
    applyAT = @(y) A' * y;
  else
    applyA = @(x) A( x, 'notransp' );
    applyAT = @(x) A( x, 'transp' );
  end

  xBar = x;
  z = applyA( x );

  if nargout > 1,  objValues = zeros( N, 1 ); end

  for optIter = 1 : N
    tmp = z + sigma * applyA( xBar );
    z = proxgConj( tmp, sigma );

    lastX = x;
    tmp = x - tau * applyAT( z );
    x = proxf( tmp );

    if nargout > 1, objValues( optIter ) = f(x) + g( applyA(x) ); end

    xBar = x + theta * ( x - lastX );
  end

  xStar = x;
end
