
function [xStar,objectiveValues,relDiffs] = gradDescent( x, gGrad, varargin )
  % [xStar,objectiveValues,relDiffs] = gradDescent( x, gGrad [, ...
  %   't', t, 'tol', tol, 'g', g, 'N', N, 'verbose', verbose ] )
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
  % t - step size (default is 1)
  % g - a function handle representing the g function; accepts a vector x
  %   as input and returns a scalar.  This is needed to calculate the
  %   objective values.
  % N - the maximum number of iterations that gradDescent will perform (default is 100)
  % verbose - if set then prints iteration information
  %
  % Outputs:
  % xStar - the optimal point
  %
  % Written by Nicholas Dwork - Copyright 2018
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp( 'Usage:  [xStar,objectiveValues,relDiffs] = gradDescent( x, gGrad [, ... ' );
    disp( '  ''t'', t, ''tol'', tol, ''g'', g, ''N'', N, ''verbose'', verbose ] ) ' );
    if nargout > 0, xStar = []; end
    if nargout > 1, objectiveValues = []; end
    return
  end

  p = inputParser;
  p.addParameter( 'g', [] );
  p.addParameter( 'N', 100, @isnumeric );
  p.addParameter( 't', 1, @isnumeric );
  p.addParameter( 'tol', [], @isnumeric );
  p.addParameter( 'verbose', 0, @isnumeric );
  p.parse( varargin{:} );
  g = p.Results.g;
  N = p.Results.N;
  t = p.Results.t;
  tol = p.Results.tol;
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

  for k = 0 : N-1
    if verbose, disp([ 'gradDescent Iteration: ', num2str(k) ]); end
    if calculateObjectiveValues == true, objectiveValues(k+1) = g(x); end

    if numel( tol ) > 0  ||  calculateObjectiveValues == true, gx = g( x ); end

    if calculateObjectiveValues == true, objectiveValues( k+1 ) = gx; end

    if numel( tol ) > 0 && tol ~= 0, lastX = x; end

    if calculateRelDiffs == true, lastX = x; end
    
    x = x - t * gGrad( x );

    if calculateRelDiffs == true
      relDiff = norm( x(:) - lastX(:) ) / norm( x(:) );
    end

    if nargout > 2, relDiffs( k+1 ) = relDiff; end

    if numel( tol ) > 0 && tol ~= 0      
      if relDiff < tol, break; end
    end
  end

  xStar = x;
end

