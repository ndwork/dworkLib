
function [xStar,objValues] = admm( x, proxf, proxg, t, varargin )
  % [xStar,objValues] = admm( x, proxf, proxg, 'A', A, 'f', f, 'g', g )
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
  p.addParameter( 'solveSys', [] );
  p.parse( varargin{:} );
  A = p.Results.A;
  f = p.Results.f;
  g = p.Results.g;
  N = p.Results.N;
  solveSys = p.Results.solveSys;

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

  z1 = x;
  z2 = applyA( x );

  if nargout > 1, objValues = numel( N, 1 ); end

  for optIter = 1 : N
    x1 = proxf( z1, t );
    x2 = proxg( z2, t );

    if nargout > 1
      if numel( f ) ==0 || numel( g ) == 0, error( 'Must supply f and g' ); end
      objValues( optIter ) = f( x ) + g( applyA( x ) );
    end

    tmp = 2*x1 - z1 + applyAT( 2*x2 - z2 );
    if numel( A ) == 0
      y1 = 0.5 * tmp;
    elseif isnumeric( A )
      ATApI = add2Diag( A' * A, 1 );
      y1 = ATApI \ tmp;
    else
      if numel( solveSys ) == 0, error( 'solveSys is required' ); end
      y1 = solveSys( tmp, applyA, applyAT );
    end
    y2 = applyA( y1 );

    z1 = z1 + y1 - x1;
    z2 = z2 + y2 - x2;
  end

  xStar = x1;
end

