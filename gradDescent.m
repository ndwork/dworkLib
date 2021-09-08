
function [xStar,objectiveValues,relDiffs] = gradDescent( x, gGrad, varargin )
  % [xStar,objectiveValues,relDiffs] = gradDescent( x, gGrad [, 't', t, 'tol', tol, ...
  %   'useMomentum', useMomentum, 'g', g, 'N', N, 'verbose', verbose ] )
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
  % g - a function handle representing the g function; accepts a vector x
  %   as input and returns a scalar.  This is needed to calculate the
  %   objective values.
  % N - the maximum number of iterations that gradDescent will perform (default is 100)
  % postOp - a function handle to a function that's called after the gradient step
  %   As an example, this could be used to implement the gradient projection method, where
  %   postOp is a handle to a function that does the projection.
  % t - step size (default is 1)
  %   Note, if g is Lipschitz with constant L, gradDescent converges with t = 1 / L.
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
    disp( '  ''t'', t, ''useMomentum'', useMomentum, ''tol'', tol, ''g'', g, ' );
    disp( '  ''N'', N, ''verbose'', verbose ] ) ' );
    if nargout > 0, xStar = []; end
    if nargout > 1, objectiveValues = []; end
    if nargout > 1, relDiffs = []; end
    return
  end

  p = inputParser;
  p.addParameter( 'g', [] );
  p.addParameter( 'N', 100, @isnumeric );
  p.addParameter( 't', 1, @isnumeric );
  p.addParameter( 'tol', [], @isnumeric );
  p.addParameter( 'useMomentum', false, @islogical );
  p.addParameter( 'verbose', 0, @islogical );
  p.parse( varargin{:} );
  g = p.Results.g;
  N = p.Results.N;
  t = p.Results.t;
  tol = p.Results.tol;
  useMomentum = p.Results.useMomentum;
  verbose = p.Results.verbose;

  if t <= 0, error('gradDescent: t0 must be greater than 0'); end

  calculateObjectiveValues = false;
  if nargout > 1
    if numel(g) == 0
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
    if calculateObjectiveValues == true, objectiveValues(k+1) = g(x); end
    
    if verbose == true
      if calculateObjectiveValues == true
        disp([ 'gradDescent Iteration: ', num2str(k), ' with objective value ', ...
          num2str( objectiveValues(k+1) ) ]);
      else
        disp([ 'gradDescent Iteration: ', num2str(k) ]);
      end
    end

    if numel( tol ) > 0  ||  calculateObjectiveValues == true, gx = g( x ); end

    if calculateObjectiveValues == true, objectiveValues( k+1 ) = gx; end

    if numel( tol ) > 0 && tol ~= 0, lastX = x; end

    if calculateRelDiffs == true, lastX = x; end

    gGradX = gGrad( x );
    
    if useMomentum == true  % Nesterov's Momentum
      lastZ = z;
      z = x - t * gGradX;
      lastMu = mu;
      mu = 0.5 * ( 1 + sqrt( 1 + 4 * mu*mu ) );
      x = z + ( lastMu / mu ) * ( z - lastZ );

    else
      x = x - t * gGradX;
    end

    if numel( postOp ) > 0
      x = postOp( x );
    end

    if calculateRelDiffs == true
      relDiff = norm( x(:) - lastX(:) ) / norm( x(:) );
    end

    if nargout > 2, relDiffs( k+1 ) = relDiff; end

    if numel( tol ) > 0 && tol ~= 0
      if relDiff < tol, break; end
    end
  end

  if nargout > 1, objectiveValues = objectiveValues( 1 : k ); end
  if nargout > 2, relDiffs = relDiffs( 1 : k ); end

  xStar = x;
end

