
function [xStar,objValues] = pdhg( x, proxf, proxgConj, tau, varargin )
  % [xStar,objValues] = pdhg( x, proxf, proxgConj, tau [, ...
  %   'A', A, 'f', f, 'g', g, 'N', N, 'normA', normA, 'sigma', sigma, ...
  %    'lambda', lambda, 'theta', theta, 'verbose', verbose, 'z', z ] )
  %
  % Implements the Primal Dual Hybrid Gradient (Chambolle-Pock) method that
  % solves problems of the form:  minimize f( x ) + g( A x )
  %
  % Inputs:
  % x - initial guess
  % proxf - a function handle for the proximal operator of
  % proxgConj - a function handle for the proximal operator of the conjugate of g
  %
  % Optional Inputs:
  % A - if A is not provided, it is assumed to be the identity
  % f - to determine the objective values, f must be provided
  % g - to determine the objective values, g must be provided
  % N - the number of iterations that ADMM will perform (default is 100)
  % normA - the matrix induced 2-norm of A.  Could be determined with norm or, for
  %   large sparse matries, estimated with normest or powerIteration.
  % lambda - relaxation parameter
  % theta - acceleration parameter
  % verbose - true or false
  % z - initial value of dual variable
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
  p.addParameter( 'normA', [], @ispositive );
  p.addParameter( 'sigma', [], @ispositive );
  p.addParameter( 'theta', 1, @(x) x >= 0 && x <= 1 );
  p.addParameter( 'lambda', 1, @(x) x >= 0 && x <= 2 );
  p.addParameter( 'verbose', false, @(x) islogical(x) || x==1 || x==0 );
  p.addParameter( 'printEvery', 1, @ispositive );
  p.addParameter( 'z', [], @isnumeric );
  p.parse( varargin{:} );
  A = p.Results.A;
  f = p.Results.f;
  g = p.Results.g;
  N = p.Results.N;
  normA = p.Results.normA;
  sigma = p.Results.sigma;
  theta = p.Results.theta;
  lambda = p.Results.lambda;
  printEvery = p.Results.printEvery;
  verbose = p.Results.verbose;
  z = p.Results.z;

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

  if numel( sigma ) == 0  && numel( A ) > 0
    if numel( normA ) == 0
      error( 'If an A is supplied, you must supply sigma or normA' );
    end
    sigma = ( 0.99 / normA^2 ) / tau;
  end

  if nargout > 1,  objValues = zeros( N, 1 ); end

  for optIter = 1 : N
    if nargout > 1
      objValues( optIter ) = f( x ) + g( applyA( x ) );
    end

    if verbose == true
      if mod(optIter,printEvery) == 0 || optIter == 1
        if nargout > 1
          disp([ 'pdhg: working on ', indx2str(optIter,N), ' of ', num2str(N), ',  ', ...
            'objective value: ', num2str( objValues( optIter ) ) ]);
        else
          disp([ 'pdhg: working on ', indx2str(optIter,N), ' of ', num2str(N) ]);
        end
      end
    end

    if optIter == 1 && numel( z ) == 0
      % z == 0 => ATz=0 => tmp == x
      z = 0;
      tmp = x;
    else
      ATz = applyAT( z );
      tmp = x - tau * ATz;
    end
    xBar = proxf( tmp, tau );

    xTmp = xBar + theta * ( xBar - x );
    AxTmp = applyA( xTmp );
    tmp = z + sigma * AxTmp;
    zBar = proxgConj( tmp, sigma );

    x = x + lambda * ( xBar - x );
    z = z + lambda * ( zBar - z );
  end

  xStar = x;
end
