
function [xStar,objectiveValues,relDiffs] = fista_wRestart( x, gGrad, proxth, varargin )
  % [xStar,objectiveValues,relDiffs] = fista_wRestart( x, g, gGrad, proxth [, ...
  %   'g', g, 'h', h, 'N', N, 'tol', tol, 'verbose', verbose ] )
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
  %
  % Optional Inputs:
  % g - a function handle representing the g function; accepts a vector x
  %     as input and returns a scalar.
  % h - a handle to the h function.  This is needed to calculate the
  %     objective values.
  % N - the number of iterations that FISTA will perform
  % t - step size (default is 1)
  % tol - if relDiff is less than this value then optimization ends
  % verbose - if set then prints fista iteration
  %
  % Outputs:
  % xStar - the optimal point
  %
  % Optional Outputs:
  % objectiveValues - value for each iteration
  % relDiffs - value for each iteration
  %
  % Written by Nicholas Dwork - Copyright 2017
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp( 'Usage:  [xStar,objectiveValues,relDiffs] = fista_wRestart( x, g, gGrad, proxth [, ... ' );
    disp( '  ''g'', g, ''h'', h, ''N'', N, ''tol'', tol, ''verbose'', verbose ] ) ' );
    if nargout > 0, xStar = [];  end
    if nargout > 1, objectiveValues = []; end
    if nargout > 2, relDiffs = []; end
    return
  end


  defaultTol = 1d-4;
  
  p = inputParser;
  p.addParameter( 'g', [] );
  p.addParameter( 'h', [] );
  p.addParameter( 'N', 100, @isnumeric );
  p.addParameter( 't', 1, @isnumeric );
  p.addParameter( 'tol', defaultTol, @(x) isnumeric(x) || numel(x) == 0 );
  p.addParameter( 'verbose', 0, @isnumeric );
  p.parse( varargin{:} );
  g = p.Results.g;
  h = p.Results.h;
  N = p.Results.N;  % total number of iterations
  t = p.Results.t;  % t0 must be greater than 0
  tol = p.Results.tol;
  verbose = p.Results.verbose;

  if t <= 0, error('fista: t0 must be greater than 0'); end
  
  calculateObjectiveValues = false;
  if nargout > 1
    if numel(h) == 0
      warning('fista.m - Cannot calculate objective values without h function handle');
    else
      objectiveValues = zeros( N, 1 );
      calculateObjectiveValues = true;
    end
  end

  calculateRelDiffs = false;
  if numel( tol ) > 0  &&  tol > 0  &&  tol < Inf
    calculateRelDiffs = true;
    relDiffs = zeros( N, 1 );
  end
  
  z = x;
  y = [];
  restarted = 0;

  k = 0;
  iter = 0;
  nRestarts = 0;
  while iter < N

    if verbose, disp([ 'FISTA wRestarting Iteration: ', num2str(iter) ]); end
    if calculateObjectiveValues > 0, objectiveValues(iter+1) = g(x) + h(x); end

    gGradZ = gGrad( z );
    x = z - t * gGradZ;

    lastY = y;
    y = proxth( x, t );
    if numel( lastY ) == 0, lastY = y; end
    
    lastZ = z;
    z = y + (k/(k+3)) * (y-lastY);

    if calculateRelDiffs == true
      relDiff = norm( z(:) - lastZ(:) ) / norm( lastZ(:) );
      relDiffs( iter + 1 ) = relDiff;
      if relDiff < tol
        break;
      end
    end

    traj = z - lastZ;
    if dotP( traj, gGradZ ) > 0 && restarted == 0
      % Restarts when trajectory and -gradient form oblique angles
      restarted = 1;
      y = [];
      k = 0;
      nRestarts = nRestarts + 1;
    else
      restarted = 0;
      k = k + 1;
    end
    iter = iter + 1;

  end

  if calculateObjectiveValues == true
    objectiveValues = objectiveValues( 1 : iter );
  end
  if calculateRelDiffs == true
    relDiffs = relDiffs( 1 : iter );
  end

  xStar = y;
  if verbose == true
    disp([ 'Fista w Restart Number of restarts: ', num2str(nRestarts) ]);
  end
end

