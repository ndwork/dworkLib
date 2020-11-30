
function [xStar,objValues] = pdhgAdaptive( x, proxf, proxgConj, tau, varargin )
  % [xStar,objValues] = pdhgAdaptive( x, proxf, proxgConj, tau [, ...
  %   'alpha0', alpha0, 'eta', eta, 'Delta', Delta', 's', s, ...  % adaptive variables
  %   'beta', beta', 'gamma', gamma, % backtracking parameters
  %   'A', A, 'f', f, 'g', g, 'innerProd', innerProd, 'N', N, 'normA', normA, ...
  %   'sigma', sigma, 'y', y, 'verbose', verbose ] )
  %
  % Implements the Adaptive Primal Dual Hybrid Gradient (Chambolle-Pock) method
  % of "Adaptive Primal-Dual Hybrid Gradient Methods for Saddle-Point Problems"
  % by Tom Goldstein et al., 2015
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
  % sigma - the second step size
  % y - initial values for dual variable
  % verbose - true or false
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
  p.addParameter( 'alpha', 0.5, @(x) x>0 && x<1 );
  p.addParameter( 'beta', 0.95, @(x) x>0 && x<1 );
  p.addParameter( 'Delta', 1.5, @(x) x>1 );
  p.addParameter( 'eta', 0.999, @(x) x>0 && x<1 );
  p.addParameter( 'f', [] );
  p.addParameter( 'g', [] );
  p.addParameter( 'gamma', 0.75, @(x) x>0 && x<1 );
  p.addParameter( 'innerProd', [] );
  p.addParameter( 'N', 100, @ispositive );
  p.addParameter( 'normA', [], @ispositive );
  p.addParameter( 'printEvery', 1, @ispositive );
  p.addParameter( 's', 1, @ispositive );
  p.addParameter( 'sigma', [], @ispositive );
  p.addParameter( 'verbose', false, @(x) islogical(x) || x==1 || x==0 );
  p.addParameter( 'y', [], @isnumeric );
  p.parse( varargin{:} );
  A = p.Results.A;
  beta = p.Results.beta;
  alpha = p.Results.alpha;
  Delta = p.Results.Delta;
  eta = p.Results.eta;
  f = p.Results.f;
  g = p.Results.g;
  gamma = p.Results.gamma;
  innerProd = p.Results.innerProd;
  N = p.Results.N;
  normA = p.Results.normA;
  printEvery = p.Results.printEvery;
  s = p.Results.s;
  sigma = p.Results.sigma;
  verbose = p.Results.verbose;
  y = p.Results.y;

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

  if numel( innerProd ) == 0
    innerProd = @(x,y) real( dotP( x, y ) );
  end

  if numel( sigma ) == 0  && numel( A ) > 0
    if numel( normA ) == 0
      error( 'If an A is supplied, you must supply sigma or normA' );
    end
    sigma = ( 0.95 / normA^2 ) / tau;
  end

  Ax = applyA( x );

  if numel( y ) == 0
    y = zeros( size( Ax ) );
  end

  if nargout > 1,  objValues = zeros( N, 1 ); end

  Ax = applyA( x );
  ATy = applyAT( y );

  for optIter = 1 : N   
    if nargout > 1
      objValues( optIter ) = f( x ) + g( Ax );
    end
    if verbose == true
      if mod(optIter,printEvery) == 0 || optIter == 1
        if nargout > 1
          disp([ 'pdhgAdaptive: working on ', indx2str(optIter,N), ' of ', num2str(N), ',  ', ...
            'objective value: ', num2str( objValues( optIter ), '%15.13f' ) ]);
        else
          disp([ 'pdhgAdaptive: working on ', indx2str(optIter,N), ' of ', num2str(N) ]);
        end
      end
    end

    lastX = x;
    lastY = y;
    lastAx = Ax;
    lastATy = ATy;

    
    % Backtracking
    while true
      tmpf = lastX - tau * lastATy;
      x = proxf( tmpf, tau );
      Ax = applyA( x );

      tmpg = lastY + sigma * ( 2 * Ax - lastAx ) ;
      y = proxgConj( tmpg, sigma );

      if optIter == 1, break; end

      bNum = 2 * tau * sigma * innerProd( Ax - lastAx, y - lastY );
      normDiffXSq = innerProd( x(:) - lastX(:), x(:) - lastX(:) );
      normDiffYSq = innerProd( y(:) - lastY(:), y(:) - lastY(:) );
      bDen = gamma * ( sigma * normDiffXSq + tau * normDiffYSq );
      b = bNum / bDen;

      if b <= 1 || bDen == 0, break; end

      % backtracking property does not hold
      tau = beta * tau / b;
      sigma = beta * sigma / b;
      
      if verbose == true
        disp(['   Backtracked New Steps tau / sigma: ', num2str(tau), ' / ', num2str(sigma) ]);
      end
    end

    ATy = applyAT( y );
    if optIter == 1, continue; end

    px = ( lastX - x ) / tau;
    py = ( lastATy - ATy );
    p = px - py;
    pRes = norm( p(:) );  % primal residual

    dy = ( lastY - y ) / sigma;
    dx = ( lastAx - Ax );
    d = dy - dx;
    dRes = norm( d(:) );  % dual residual

    if pRes > s * dRes * Delta
      % primal residual is large
      tau = tau / ( 1 - alpha );
      sigma = sigma * ( 1 - alpha );
      alpha = alpha * eta;
      if verbose == true
        disp(['   New Steps tau / sigma: ', num2str(tau), ' / ', num2str(sigma) ]);
      end

    elseif pRes < s * dRes / Delta
      % dual residual is large
      tau = tau * ( 1 - alpha );
      sigma = sigma / ( 1 - alpha );
      alpha = alpha * eta;
      if verbose == true
        disp(['   New Steps tau / sigma: ', num2str(tau), ' / ', num2str(sigma) ]);
      end

    else
      % residuals are similar
      % do not change step sizes
    end
  end

  xStar = x;
end
