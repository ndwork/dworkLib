
function [xStar,objectiveValues] = proxGrad( x, g, gGrad, proxth, varargin )
  % [xStar,objectiveValues] = proxGrad( x, g, gGrad, proxth [, ...
  %   't', 'h', h, 'N', N, 'verbose', verbose ] )
  %
  % This function implements the proximal gradient method.
  % This method finds the x that minimizes functions of form g(x) + h(x) where
  % g is differentiable and h has a simple proximal operator.
  %
  % Inputs:
  % x - the starting point
  % g - a function handle representing the g function; accepts a vector x
  %     as input and returns a scalar.
  % gGrad - a function handle representing the gradient function of g;
  %     input: the point to evaluation, output: the gradient vector
  % proxth - the proximal operator of the h function (with parameter t);
  %     two inputs: the vector and the scalar value of the parameter t
  %
  % Optional Inputs:
  % t - step size (default is 1)
  % h - a handle to the h function.  This is needed to calculate the
  %     objective values.
  % N - the number of iterations that proxGrad will perform
  % verbose - if set then prints iteration information
  %
  % Outputs:
  % xStar - the optimal point
  %
  % Written by Nicholas Dwork - Copyright 2018
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addParameter( 't', 1, @isnumeric );
  p.addParameter( 'h', [] );
  p.addParameter( 'N', 100, @isnumeric );
  p.addParameter( 'verbose', 0, @isnumeric );
  p.parse( varargin{:} );
  h = p.Results.h;
  N = p.Results.N;  % total number of iterations
  t = p.Results.t;  % t0 must be greater than 0
  verbose = p.Results.verbose;

  if t <= 0, error('fista: t0 must be greater than 0'); end
  
  calculateObjectiveValues = 0;
  if nargout > 1
    if numel(h) == 0
      warning('fista.m - Cannot calculate objective values without h function handle');
    else
      objectiveValues = zeros(N,1);
      calculateObjectiveValues = 1;
    end
  end

  for k=0:N-1
    if verbose, disp([ 'proxGrad Iteration: ', num2str(k) ]); end
    if numel(calculateObjectiveValues) > 0, objectiveValues(k+1) = g(x) + h(x); end

    y = x - t * gGrad( x );
    x = proxth( y, t );
  end

  xStar = x;
end

