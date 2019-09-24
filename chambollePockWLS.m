
function [xStar,objValues] = chambollePockWLS( x, proxf, proxgConj, varargin )
  % [xStar,objValues] = chambollePockWLS( x, proxf, proxgConj [, ...
  %   'N', N, 'A', A, 'f', f, 'g', g, 'mu', mu, 'tau', tau, 'theta', theta, ...
  %   'verbose', verbose ] )
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
  p.addParameter( 'beta', 1, @ispositive );
  p.addParameter( 'delta', 0.99, @(x) x>0 && x<1 );
  p.addParameter( 'doCheckAdjoint', false, @islogical );
  p.addParameter( 'f', [] );
  p.addParameter( 'g', [] );
  p.addParameter( 'mu', 0.7, @(x) x>0 && x<1 );
  p.addParameter( 'N', 100, @ispositive );
  p.addParameter( 'tau', 1, @ispositive );
  p.addParameter( 'theta', 1, @ispositive );
  p.addParameter( 'verbose', false, @islogical );
  p.parse( varargin{:} );
  A = p.Results.A;
  beta = p.Results.beta;
  delta = p.Results.delta;
  doCheckAdjoint = p.Results.doCheckAdjoint;
  f = p.Results.f;
  g = p.Results.g;
  mu = p.Results.mu;
  N = p.Results.N;
  tau = p.Results.tau;
  theta = p.Results.theta;
  verbose = p.Results.verbose;

  if numel( A ) == 0
    applyA = @(x) x;
    applyAT = @(x) x;
    y = x;
  elseif isnumeric( A )
    applyA = @(x) A * x;
    applyAT = @(y) A' * y;
    y = A * x;
  else
    applyA = @(x) A( x, 'notransp' );
    applyAT = @(x) A( x, 'transp' );
    y = applyA( x );
  end

  if doCheckAdjoint == true
    [adjointCheckPassed,adjCheckErr] = checkAdjoint( x, applyA, applyAT );
    if ~adjointCheckPassed, error([ 'checkAdjoint failed with error ', num2str(adjCheckErr) ]); end
  end

  if nargout > 1,  objValues = zeros( N, 1 ); end

  for optIter = 1 : N
    if verbose == true
      disp([ 'chambollePockWLS: working on ', num2str(optIter), ' of ', num2str(N) ]);
    end

    lastX = x;
    tmp = lastX - tau * applyAT( y );
    x = proxf( tmp, tau );

    if nargout > 1
      objValues( optIter ) = f( x ) + g( applyA( x ) );
    end

    lastTau = tau;
    tau = tau * sqrt( 1 + theta );

    diffx = x - lastX;

    lastY = y;
    while true
      theta = tau / lastTau;

      xBar = x + theta * ( diffx );

      betaTau = beta * tau;
      tmp = lastY + betaTau * applyA( xBar );
      y = proxgConj( tmp, betaTau );

      diffy = y - lastY;
      ATdiffy = applyAT( diffy );

      if tau * sqrt( beta ) * norm( ATdiffy(:) ) <= delta * norm( diffy(:) )
        break
      end
      tau = mu * tau;
    end
  end

  xStar = x;
end
