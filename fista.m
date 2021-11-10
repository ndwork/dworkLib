
function [xStar,objectiveValues,relDiffs] = fista( x, gGrad, proxth, varargin )
  % [xStar,objectiveValues,relDiffs] = fista( x, gGrad, proxth [, ...
  %   'g', g, 'h', h, 'N', N, 't', t, 'verbose', verbose ] )
  %
  % This function implements the FISTA optimization algorithm
  % FISTA finds the x that minimizes functions of form g(x) + h(x) where
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
  % h - a handle to the h function.  This is needed to calculate the
  %     objective values.
  % N - the number of iterations that FISTA will perform (default is 100)
  % t - step size (default is 1)
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

  defaultN = 100;

  p = inputParser;
  p.addParameter( 'g', [] );
  p.addParameter( 'h', [] );
  p.addParameter( 'N', defaultN, @(x) ispositive(x) || numel(x)==0 );
  p.addParameter( 'printEvery', 1, @ispositive );
  p.addParameter( 't', 1, @isnumeric );
  p.addParameter( 'verbose', 0, @(x) isnumeric(x) || islogical(x) );
  p.parse( varargin{:} );
  g = p.Results.g;
  h = p.Results.h;
  N = p.Results.N;  % total number of iterations
  printEvery = p.Results.printEvery;  % display result printEvery iterations
  t = p.Results.t;  % t0 must be greater than 0
  verbose = p.Results.verbose;

  if numel( N ) == 0, N = defaultN; end

  if t <= 0, error('fista: t0 must be greater than 0'); end

  calculateObjectiveValues = 0;
  if nargout > 1
    if numel(h) == 0 || numel(g) == 0
      error( 'fista.m - Cannot calculate objective values without g and h function handles' );
    else
      objectiveValues = zeros( N, 1 );
      calculateObjectiveValues = 1;
    end
  end

  calculateRelDiffs = false;
  if nargout > 2 || verbose == true
    relDiffs = zeros( N, 1 );
    calculateRelDiffs = true;
  end
  
  z = x;
  y = [];

  for k = 0 : N-1
    x = z - t * gGrad( z );

    lastY = y;
    y = proxth( x, t );
    if numel( y ) == 0, y = lastY; end

    if calculateObjectiveValues > 0, objectiveValues(k+1) = g(y) + h(y); end

    if calculateRelDiffs == true
      relDiff = norm( z(:) - lastZ(:) ) / norm( lastZ(:) );
      relDiffs(k+1) = relDiff;
    end

    if verbose>0 && mod( k, printEvery ) == 0
      formatString = ['%', num2str(ceil(log10(N))), '.', num2str(ceil(log10(N))), 'i' ];
      verboseString = [ 'FISTA Iteration: ', num2str(k,formatString) ];
      if calculateObjectiveValues > 0
        verboseString = [ verboseString, ',  objective: ', num2str( objectiveValues(k+1) ) ];   %#ok<AGROW>
      end
      if calculateRelDiffs > 0
        verboseString = [ verboseString, ',  relDiff: ', num2str( relDiffs(k+1) ) ];   %#ok<AGROW>
      end
      disp( verboseString );
    end

    z = y + ( k / (k+3) ) * ( y - lastY );
  end

  xStar = y;
end

