
function [xStar,objectiveValues,relDiffs, restarts] = fista_wRestart( x, gGrad, proxth, q, varargin )
  % [xStar,objectiveValues,relDiffs] = fista_wRestart( x, gGrad, proxth [, ...
  %   'g', g, 'h', h, 'N', N, 't', t, 'tol', tol, 'verbose', verbose ] )
  %
  % This function implements the FISTA optimization algorithm
  % FISTA finds the x that minimizes functions of form g(x) + h(x) where
  % g is differentiable and h has a simple proximal operator.
  % Written according to "Adaptive Restart for Accelerated Gradient Schemes"
  % by Donoghue and Candes
  %
  % Inputs:
  % x - the starting point
  % gGrad - a function handle representing the gradient function of g;
  %     input: the point to evaluation, output: the gradient vector
  % proxth - the proximal operator of the h function (with parameter t);
  %     two inputs: the vector and the scalar value of the parameter t
  % q - an integer specifying how many iterations before restart
  %
  % Optional Inputs:
  % g - a function handle representing the g function; accepts a vector x
  %     as input and returns a scalar.
  % h - a handle to the h function.  This is needed to calculate the
  %     objective values.
  % N - the number of iterations that FISTA will perform (default is 100)
  % t - step size (default is 1)
  % tol - runs until relDiff is less than tol or until iterations is N
  % verbose - if set then prints fista iteration
  %
  % Outputs:
  % xStar - the optimal point
  %
  % Written by Nicholas Dwork - Copyright 2017
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp( 'Usage:  [xStar,objectiveValues,relDiffs] = fista_wRestart( x, gGrad, proxth, q [, ... ' );
    disp( '  ''g'', g, ''h'', h, ''N'', N, ''tol'', tol, ''verbose'', verbose ] ) ' );
    if nargout > 0, xStar = [];  end
    if nargout > 1, objectiveValues = []; end
    if nargout > 2, relDiffs = []; end
    return
  end

  defaultN = 100;

  p = inputParser;
  p.addParameter( 'g', [] );
  p.addParameter( 'h', [] );
  p.addParameter( 'N', defaultN, @(x) ispositive(x) || numel(x)==0 );
  p.addParameter( 'printEvery', 1, @ispositive );
  p.addParameter( 't', 1, @isnumeric );
  p.addParameter( 'tol', 1d-4, @(x) numel(x) == 0  ||  x >= 0 );
  p.addParameter( 'verbose', 0, @(x) isnumeric(x) || islogical(x) );
  p.parse( varargin{:} );
  g = p.Results.g;
  h = p.Results.h;
  N = p.Results.N;  % total number of iterations
  printEvery = p.Results.printEvery;  % display result printEvery iterations
  t = p.Results.t;  % t0 must be greater than 0
  tol = p.Results.tol;
  verbose = p.Results.verbose;

  if numel( N ) == 0, N = defaultN; end

  if t <= 0, error('fista_wRestart: t0 must be greater than 0'); end

  if nargout > 1, objectiveValues = zeros( N, 1 ); end
  if nargout > 2, relDiffs = zeros( N, 1 ); end
  if nargout > 3, restarts = zeros( N, 1 ); end

  calculateObjectiveValues = 0;
  if nargout > 1
    if numel(h) == 0 || numel(g) == 0
      error( 'fista_wRestart.m - Cannot calculate objective values without g and h function handles' );
    else
      calculateObjectiveValues = true;
    end
  end

  calculateRelDiffs = false;
  if nargout > 2  ||  verbose == true  ||  ( numel(tol) > 0 && tol < Inf )
    calculateRelDiffs = true;
  end

  z = x;
  y = x;
  maxAbsZ = max( abs( z(:) ) );

  k = -1;
  iter = 0;
  while iter < N
    k = k + 1;
    iter = iter + 1;

    x = z - t * gGrad( z );

    lastY = y;
    y = proxth( x, t );

    restarted = false;
    if mod( iter, q ) == 0
      k = 0;  % Restart by eliminating momentum
      restarted = true;
    end

    lastZ = z;
    lastMaxAbsZ = maxAbsZ;
    z = y + ( k / (k+3) ) * ( y - lastY );

    maxAbsZ = max( abs( z(:) ) );

    if calculateObjectiveValues > 0, objectiveValues(iter) = g(z) + h(z); end

    if calculateRelDiffs == true
      relDiff = norm( z(:) - lastZ(:) ) / norm( lastZ(:) );
    end
    if nargout > 2  &&  iter > 1, relDiffs(iter) = relDiff; end

    if verbose>0 && mod( iter, printEvery ) == 0
      formatString = ['%', num2str(ceil(log10(N))), '.', num2str(ceil(log10(N))), 'i' ];
      verboseString = [ 'FISTA wRestart Iteration: ', num2str(iter,formatString) ];
      if calculateObjectiveValues > 0
        verboseString = [ verboseString, ',  objective: ', num2str( objectiveValues(iter) ) ];   %#ok<AGROW>
      end
      if calculateRelDiffs > 0
        verboseString = [ verboseString, ',  relDiff: ', num2str( relDiff ) ];   %#ok<AGROW>
      end
      if restarted == true
        verboseString = [ verboseString, ' - RESTARTED' ];   %#ok<AGROW>
      end
      disp( verboseString );
    end

    if numel(tol) > 0  &&  tol < Inf  &&  relDiff < tol
      break;
    end
    if maxAbsZ == 0  && lastMaxAbsZ == 0
      break;
    end
  end

  if nargout > 1, objectiveValues = objectiveValues( 1 : iter ); end
  if nargout > 2, relDiffs = relDiffs( 1 : iter ); end
  if nargout > 3, restarts = restarts( 1 : iter ); end

  xStar = y;
end

