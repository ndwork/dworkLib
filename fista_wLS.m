
function [xStar,objectiveValues] = fista_wLS( x, g, gGrad, proxth, varargin )
  % [xStar,optValue] = fista( x, g, gGrad, proxth [, ...
  %   'h', h, 'N', N, 'r', r, 's', s, 't0', t0, 'verbose', verbose ] )
  %
  % This function implements the FISTA optimization algorithm with line
  % search as described in "Fraction-variant beam orientation optimization
  % for non-coplanar IMRT" by O'Connor et al. (2017)
  %
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
  % r - the backtracking line search parameter; must be between 0 and 1
  %     (default is 0.5)
  % s - the scaling parameter each FISTA iteration; must be greater than 1
  %     (default is 1.25)
  % t0 - initial step size (default is 1)
  % verbose - if set then prints fista iteration
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

  defaultN = 100;

  p = inputParser;
  p.addParameter( 'h', [] );
  p.addParameter( 'N', defaultN, @(x) ispositive(x) || numel(x) == 0 );
  p.addParameter( 'printEvery', 1, @ispositive );
  p.addParameter( 'r', 0.5, @ispositive );
  p.addParameter( 's', 1.25, @ispositive );
  p.addParameter( 't0', 1, @ispositive );
  p.addParameter( 'verbose', false, @(x) isnumeric(x) || islogical(x) );
  p.parse( varargin{:} );
  h = p.Results.h;
  N = p.Results.N;  % total number of iterations
  printEvery = p.Results.printEvery;  % display result printEvery iterations
  r = p.Results.r;  % r must be between 0 and 1
  s = p.Results.s;  % s must be greater than 1
  t0 = p.Results.t0;  % t0 must be greater than 0
  verbose = p.Results.verbose;

  if numel( N ) == 0, N = defaultN; end

  if r <= 0 || r >= 1, error('fista: r must be in (0,1)'); end
  if s <= 1, error('fista: s must be greater than 1'); end
  if t0 <= 0, error('fista: t0 must be greater than 0'); end
  
  calculateObjectiveValues = 0;
  if nargout > 1
    if numel(h) == 0
      warning('fista.m - Cannot calculate objective values without h function handle');
    else
      objectiveValues = zeros(N,1);
      calculateObjectiveValues = 1;
    end
  end


  t = t0;
  v = x;
  theta = 1;

  for k=0:N-1
    lastX = x;
    lastT = t;
    t = s * lastT;
    lastTheta = theta;

    while true
      if k==0
        theta = 1;
      else
        a = lastT;  b = t*lastTheta*lastTheta;  c = -b;
        theta = ( -b + sqrt( b*b - 4*a*c ) ) / ( 2*a );
      end
      y = (1-theta) * lastX + theta * v;

      Dgy = gGrad( y );
      x = proxth( y - t * Dgy, t );

      breakThresh = g(y) + dotP(Dgy,x-y) + (1/(2*t)) * norm( x(:) - y(:), 2 )^2;
      if g(x) <= breakThresh, break; end
      t = r*t;
      if verbose>1 && mod( k, printEvery ) == 0
        disp([ '  Step size change to: ', num2str(t) ]);
      end
    end

    if calculateObjectiveValues > 0, objectiveValues(k+1) = g(x) + h(x); end
    if verbose>0 && mod( k, printEvery ) == 0
      formatString = ['%', num2str(ceil(log10(N))), '.', num2str(ceil(log10(N))), 'i' ];
      verboseString = [ 'FISTA (with Line Search) Iteration: ', num2str(k,formatString) ];
      if calculateObjectiveValues > 0
        verboseString = [ verboseString, ',  objective: ', ...
          num2str( objectiveValues(k+1) ) ];   %#ok<AGROW>
      end
      disp( verboseString );
    end

    v = x + (1/theta) * ( x - lastX );
  end

  xStar = x;
end

