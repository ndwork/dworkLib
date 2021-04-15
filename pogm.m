
function [xStar,objectiveValues] = pogm( x, gGrad, proxth, varargin )
  % [xStar,objectiveValues] = pogm( x, g, gGrad, proxth [, N, ...
  %   'g', g, 'h', h, 't', t, 'verbose', verbose ] )
  %
  % This function implements the POGM optimization algorithm
  % POGM finds the x that minimizes functions of form g(x) + h(x) where
  % g is differentiable and h has a simple proximal operator.
  % Though the same order as FISTA, it converges faster when the number of
  % iterations is fixed.
  %
  % This is written according to Fig. 3 in "Optimization Methods for MR
  % image reconstruction (Long Version)" by Jeff Fessler
  %
  % Inputs:
  % x - the starting point
  % gGrad - a function handle representing the gradient function of g;
  %     input: the point to evaluation, output: the gradient vector
  % proxth - the proximal operator of the h function (with parameter t);
  %     two inputs: the vector and the scalar value of the parameter t
  % N - the number of iterations that FISTA will perform (default is 100)
  %
  % Optional Inputs:
  % g - a function handle representing the g function; accepts a vector x
  %     as input and returns a scalar.  This is needed to calculate the
  %     objective values.
  % h - a handle to the h function.  This is needed to calculate the
  %     objective values.
  % t - step size (default is 1)
  % verbose - if set then prints pogm iteration
  %
  % Outputs:
  % xStar - the optimal point
  %
  % Written by Nicholas Dwork - Copyright 2021
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addOptional( 'N', 100, @ispositive );
  p.addParameter( 'g', [] );
  p.addParameter( 'h', [] );
  p.addParameter( 'printEvery', 1, @ispositive );
  p.addParameter( 't', 1, @isnumeric );
  p.addParameter( 'verbose', 0, @(x) isnumeric(x) || islogical(x) );
  p.parse( varargin{:} );
  g = p.Results.g;
  h = p.Results.h;
  N = p.Results.N;
  printEvery = p.Results.printEvery;  % display result printEvery iterations
  t = p.Results.t;  % t0 must be greater than 0
  verbose = p.Results.verbose;

  if t <= 0, error('fista: t0 must be greater than 0'); end

  calculateObjectiveValues = 0;
  if nargout > 1
    if numel(h) == 0
      error( 'fista.m - Cannot calculate objective values without h function handle' );
    else
      objectiveValues = zeros(N,1);
      calculateObjectiveValues = 1;
    end
  end

  theta = 1;
  w = x;
  z = x;

  for k = 0:N-1

    lastTheta = theta;
    if k < N-1
      theta = 0.5 * ( 1 + sqrt( 4 * theta * theta + 1 ) );
    else
      theta = 0.5 * ( 1 + sqrt( 8 * theta * theta + 1 ) );
    end
    gamma = t * ( 2*lastTheta + theta - 1 ) / theta;

    lastW = w;
    w = x - t * gGrad( x );

    z = w ...
      + ( lastTheta - 1 ) / theta * ( w - lastW ) ...
      + lastTheta / theta * ( w - x ) ...
      + t * ( lastTheta - 1 ) / ( gamma * theta ) * ( z - x );

    x = proxth( z, gamma );
    if calculateObjectiveValues > 0, objectiveValues(k+1) = g(x) + h(x); end

    if verbose>0 && mod( k, printEvery ) == 0
      formatString = ['%', num2str(ceil(log10(N))), '.', num2str(ceil(log10(N))), 'i' ];
      verboseString = [ 'POGM Iteration: ', num2str(k,formatString) ];
      if calculateObjectiveValues > 0
        verboseString = [ verboseString, ',  objective: ', ...
          num2str( objectiveValues(k+1), 15 ) ];   %#ok<AGROW>
      end
      disp( verboseString );
    end

  end

  xStar = x;
end
