
function [xStar,objectiveValues] = gradDescent( x, gGrad, varargin )
  % [xStar,objectiveValues] = proxGrad( x, gGrad [, ...
  %   't', 'g', g, 'N', N, 'verbose', verbose ] )
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
  % N - the number of iterations that gradDescent will perform (default is 100)
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

  p = inputParser;
  p.addParameter( 't', 1, @isnumeric );
  p.addParameter( 'N', 100, @isnumeric );
  p.addParameter( 'verbose', 0, @isnumeric );
  p.parse( varargin{:} );
  N = p.Results.N;  % total number of iterations
  t = p.Results.t;  % t0 must be greater than 0
  verbose = p.Results.verbose;

  if t <= 0, error('gradDescent: t0 must be greater than 0'); end

  calculateObjectiveValues = 0;
  if nargout > 1
    if numel(g) == 0
      warning('gradDescent.m - Cannot calculate objective values without h function handle');
    else
      objectiveValues = zeros(N,1);
      calculateObjectiveValues = 1;
    end
  end

  for k=0:N-1
    if verbose, disp([ 'gradDescent Iteration: ', num2str(k) ]); end
    if numel(calculateObjectiveValues) > 0, objectiveValues(k+1) = g(x); end

    x = x - t * gGrad( x );
  end

  xStar = x;
end

