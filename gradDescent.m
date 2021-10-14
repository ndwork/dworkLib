
function [xStar,objectiveValues,relDiffs] = gradDescent( x, gGrad, varargin )
  % [xStar,objectiveValues,relDiffs] = gradDescent( x, gGrad [, 'alpha', alpha, 'beta', beta, ...
  %   't', t, 'tol', tol, 'useLineSearch', useLineSearch, 'useMomentum', useMomentum, ...
  %   'g', g, 'N', N, 'nMaxLineSearchIter', nMaxLineSearchIter, 'verbose', verbose ] )
  %
  % This function implements the gradient descent method.
  % This method finds the x that minimizes a (sub)differentiable function g
  %
  % Inputs:
  % x - the starting point
  % gGrad - a function handle representing the (sub)gradient function of g;
  %   input: the point to evaluation, output: the gradient vector
  %
  % Optional Inputs:
  % alpha - line search parameter for 
  % beta - line search parameter for step size reduction (default is 0.8)
  % g - a function handle representing the g function; accepts a vector x
  %   as input and returns a scalar.  This is needed to calculate the
  %   objective values.
  % N - the maximum number of iterations that gradDescent will perform (default is 100)
  % nMaxLineSearchIter - maximum number of line search iterations (default is 10)
  % postOp - a function handle to a function that's called after the gradient step
  %   As an example, this could be used to implement the gradient projection method, where
  %   postOp is a handle to a function that does the projection.
  % t - step size (default is 1)
  %   Note, if g is Lipschitz with constant L, gradDescent converges with t = 1 / L.
  % tau - line search parameter for step size increase (default is 1.1)
  % useLineSearch - implements Armijo line search
  % useMomentum - logical parameter that specifies whether or not to
  %   use momentum (default is false)
  % verbose - if set then prints iteration information
  %
  % Outputs:
  % xStar - the last point in the trajectory
  %
  % Written by Nicholas Dwork - Copyright 2018
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  if nargin < 1
    disp( 'Usage:  [xStar,objectiveValues,relDiffs] = gradDescent( x, gGrad [, ... ' );
    disp( '  ''alpha'', alpha, ''beta'', beta, ''t'', t, ''useLineSearch'', useLineSearch, ... ' );
    disp( '  ''useMomentum'', useMomentum, ''tol'', tol, ''g'', g, ''N'', N, ... ' );
    disp( '  ''verbose'', verbose ] ) ' );
    if nargout > 0, xStar = []; end
    if nargout > 1, objectiveValues = []; end
    if nargout > 1, relDiffs = []; end
    return
  end

  p = inputParser;
  p.addParameter( 'alpha', 0.3, @(x) x > 0 && x < 0.5 );
  p.addParameter( 'beta', 0.8, @(x) x > 0 && x < 1 );
  p.addParameter( 'g', [] );
  p.addParameter( 'N', 100, @isnumeric );
  p.addParameter( 'nMaxLineSearchIter', 100, @ispositive );
  p.addParameter( 'postOp', [] );
  p.addParameter( 'printEvery', 1, @isnumeric );
  p.addParameter( 't', 1, @isnumeric );
  p.addParameter( 'tau', 1.1, @(x) x >= 1 );
  p.addParameter( 'tol', [], @isnumeric );
  p.addParameter( 'useLineSearch', false, @islogical );
  p.addParameter( 'useMomentum', false, @islogical );
  p.addParameter( 'verbose', 0, @islogical );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  beta = p.Results.beta;
  g = p.Results.g;
  N = p.Results.N;
  nMaxLineSearchIter = p.Results.nMaxLineSearchIter;
  postOp = p.Results.postOp;
  printEvery = p.Results.printEvery;
  t = p.Results.t;
  tau = p.Results.tau;
  tol = p.Results.tol;
  useLineSearch = p.Results.useLineSearch;
  useMomentum = p.Results.useMomentum;
  verbose = p.Results.verbose;

  if t <= 0, error('gradDescent: t0 must be greater than 0'); end

  if useLineSearch == true
    if numel( g ) == 0
      error( 'gradDescent.m - line search requires g' );
    end
    if useMomentum == true
      error( 'Line search is not currently implemented with momentum.' );
    end
  end

  calculateObjectiveValues = false;
  if nargout > 1
    if numel( g ) == 0
      error( 'gradDescent.m - Cannot calculate objective values without g function handle' );
    end
    objectiveValues = zeros( N, 1 );
    calculateObjectiveValues = true;
  end

  calculateRelDiffs = false;
  if nargout > 2
    relDiffs = zeros( N, 1 );
    calculateRelDiffs = true;
  end

  if useMomentum == true
    mu = 1;
    z = x;
  end

  for k = 0 : N-1
    if numel( tol ) > 0  ||  calculateObjectiveValues == true, gx = g( x ); end

    if calculateObjectiveValues == true, objectiveValues( k+1 ) = gx; end

    if ( numel( tol ) > 0 && tol ~= 0 ) || ( calculateRelDiffs == true ), lastX = x; end

    gGradX = gGrad( x );

    if useMomentum == true  % Nesterov's Momentum
      lastZ = z;
      z = x - t * gGradX;
      lastMu = mu;
      mu = 0.5 * ( 1 + sqrt( 1 + 4 * mu*mu ) );
      x = z + ( lastMu / mu ) * ( z - lastZ );

    else

      if useLineSearch == true
        gx = g( x );
        lineSearchIter = 0;
        while true  &&  lineSearchIter <= nMaxLineSearchIter 
          lineSearchIter = lineSearchIter + 1;
          xNew = x - t * gGradX;
          if g( xNew ) < gx - alpha * t * norm( gGradX(:) )^2
            break;
          end
          t = beta * t;
        end
        t = tau * t;
        x = xNew;
      else
        x = x - t * gGradX;
      end
    end

    if numel( postOp ) > 0
      x = postOp( x );
    end

    if calculateRelDiffs == true
      relDiff = norm( x(:) - lastX(:) ) / norm( x(:) );
    end
    if nargout > 2, relDiffs( k+1 ) = relDiff; end

    if verbose == true  &&  mod( k+1, printEvery ) == 0
      verboseStr = [ 'gradDescent Iteration: ', num2str(k) ];
      if calculateObjectiveValues == true
        verboseStr = [ verboseStr, ' with objective value ', num2str( objectiveValues(k+1) ) ];   %#ok<AGROW>
      end 
      if calculateRelDiffs == true && k > 0
        verboseStr = [ verboseStr, ' with relDiff ', num2str( relDiff ) ];   %#ok<AGROW>
      end
      disp( verboseStr );
    end

    if numel( tol ) > 0 && tol ~= 0
      if relDiff < tol, break; end
    end
  end

  if nargout > 1, objectiveValues = objectiveValues( 1 : k ); end
  if nargout > 2, relDiffs = relDiffs( 1 : k ); end

  xStar = x;
end

