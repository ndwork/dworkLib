
function [xStar,objectiveValues,relDiffs] = fista_wLS( x, g, gGrad, proxth, varargin )
  % [xStar,objValues,relDiffs] = fista_wLS( x, g, gGrad, proxth [, ...
  %   'h', h, 'innerProd', innerProd, 'minStep', minStep, 'N', N, 'r', r, 's', s, 't0', t0, ...
  %   'restart', true/false, 'tol', tol, 'verbose', verbose ] )
  %
  % This function implements the FISTA optimization algorithm with line
  % search as described in "Fraction-variant beam orientation optimization
  % for non-coplanar IMRT" by O'Connor et al. (2017)
  % Restart was implemented according to "Adaptive Restart for Accelerated
  % Gradient Schemes" by Donoghue and Candes.
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
  % h - a handle to the h function.  This is needed to calculate the objective values.
  % innerProd - function handle to the inner product.
  %     ( default is real( dotP( x, y ) ) )
  %     Note: for complex vectors, innerProd should be real( dotP( x, y ) )
  % minStep - the smallest the step size can be
  % N - the number of iterations that FISTA will perform (default is 100)
  % r - the backtracking line search parameter; must be between 0 and 1
  %     (default is 0.5)
  % s - the scaling parameter each FISTA iteration; must be greater than 1
  %     (default is 1.25)
  % t0 - initial step size (default is 1)
  % tol - if the relative error is below this tolerance then fista_wLS returns
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

  if nargin < 1
    disp( 'Usage:  ' );
    disp( '  [xStar,objValues] = fista_wLS( x, g, gGrad, proxth [, ... ');
    disp( '    ''h'', h, ''innerProd'', innerProd, ''N'', N, ''r'', r, ''s'', s, ''t0'', t0, ... ');
    disp( '    ''restart'', true/false''verbose'', verbose ] ' );
    if nargout > 0, xStar = []; end
    if nargout > 1, objectiveValues = []; end
    if nargout > 2, relDiffs = []; end
    return
  end

  defaultN = 100;
  defaultSubIterThresh = 100;

  p = inputParser;
  p.addParameter( 'gradNorm', Inf, @ispositive );
  p.addParameter( 'h', [] );
  p.addParameter( 'innerProd', [] );
  p.addParameter( 'minStep', 0, @(x) x >= 0 );
  p.addParameter( 'N', defaultN, @(x) ispositive(x) || numel(x) == 0 );
  p.addParameter( 'printEvery', 1, @ispositive );
  p.addParameter( 'r', 0.5, @ispositive );
  p.addParameter( 'restart', false, @islogical );
  p.addParameter( 's', 1.25, @ispositive );
  p.addParameter( 'subIterThresh', defaultSubIterThresh, @ispositive );
  p.addParameter( 't0', 1, @ispositive );
  p.addParameter( 'tol', 1d-8, @(x) x >= 0 );
  p.addParameter( 'verbose', false, @(x) isnumeric(x) || islogical(x) );
  p.parse( varargin{:} );
  gradNorm = p.Results.gradNorm;
  h = p.Results.h;
  innerProd = p.Results.innerProd;
  minStep = p.Results.minStep;
  N = p.Results.N;  % total number of iterations
  printEvery = p.Results.printEvery;  % display result printEvery iterations
  r = p.Results.r;  % r must be between 0 and 1
  restart = p.Results.restart;
  s = p.Results.s;  % s must be greater than 1
  subIterThresh = p.Results.subIterThresh;
  t0 = p.Results.t0;  % t0 must be greater than 0
  tol = p.Results.tol;
  verbose = p.Results.verbose;

  if numel( N ) == 0, N = defaultN; end
  if numel( innerProd ) == 0
    innerProd = @(x,y) real( dotP( x, y ) );
  end

  if r <= 0 || r >= 1, error('fista: r must be in (0,1)'); end
  if s <= 1, error('fista: s must be greater than 1'); end
  if t0 <= 0, error('fista: t0 must be greater than 0'); end

  calculateObjectiveValues = false;
  if nargout > 1
    if numel(h) == 0
      warning('fista_wLS.m - Cannot calculate objective values without h function handle');
    else
      objectiveValues = zeros( N+1, 1 );
      calculateObjectiveValues = true;
    end
  end

  calculateRelDiffs = false;
  if ( numel( tol ) ~= 0 && tol > 0 ) || nargout > 2
    calculateRelDiffs = true;
    relDiffs = zeros( N, 1 );
  end

  if calculateObjectiveValues == true
    gx = g( x );
    hx = h( x );
    objectiveValues(1) = gx + hx;
  end
  t = t0 / s;
  v = x;
  theta = 1;
  k = 0;

  for iter = 1 : N
    if verbose>0 && iter>0 && mod( iter, printEvery ) == 0
      verboseString = [ 'FISTA (with Line Search) Iteration: ', indx2str(iter,N), ...
        ' of ', num2str(N) ];
      if calculateObjectiveValues == true
        verboseString = [ verboseString, ',  objective: ', ...
          num2str( objectiveValues(iter), 15 ), '   parts: ', num2str(gx), ...
          ',  ', num2str( hx ) ];   %#ok<AGROW>
      end
      if calculateRelDiffs == true  &&  iter > 1
        verboseString = [ verboseString, ', relDiff: ', num2str( relDiffs( iter-1 ) ) ];   %#ok<AGROW>
      end
      disp( verboseString );
    end

    lastX = x;
    lastT = t;
    t = s * lastT;
    lastTheta = theta;

    subIter = 0;
    while true
      subIter = subIter + 1;
      if verbose == true && mod( iter, printEvery ) == 0
        disp([ '     sub iteration: ', num2str( subIter ), '  t: ', num2str(t) ]);
      end

      if t < minStep, t = minStep; end

      if k==0
        theta = 1;
      else
        a = lastT;  b = t*lastTheta*lastTheta;  c = -b;
        %theta = ( -b + sqrt( b*b - 4*a*c ) ) / ( 2*a );
        [root1,root2] = quadRoots( a, b, c );  % Numerically stable solution to quadratic
        theta = max( root1, root2 );
      end
      y = (1-theta) * lastX + theta * v;

      Dgy = gGrad( y );
      x = proxth( y - t * Dgy, t );
      gx = g( x );
      if calculateObjectiveValues == true
        hx = h( x );
        objectiveValues( iter + 1 ) = gx + hx;
      end

      gy = g( y );
      innerProdResult = innerProd( Dgy, x-y );
      normDiffSq = innerProd( x-y, x-y );
      breakThresh = gy + innerProdResult + (1/(2*t)) * normDiffSq;

      if ( gx <= breakThresh ) || ...
         ( subIter >= subIterThresh ) || ...
         ( t <= minStep )
       break;
      end

      if t < 1 / gradNorm, break; end

      t = r*t;
      if verbose>1 && mod( k, printEvery ) == 0
        disp([ '  Step size change to: ', num2str(t) ]);
      end
    end

    if calculateRelDiffs == true
      xNorm = sqrt( innerProd( x, x ) );
      diffNorm = sqrt( innerProd( x - lastX, x - lastX ) );
      relDiff = diffNorm / xNorm;

      if verbose == true
        disp([ '  Relative error: ', num2str( relDiff ) ]);
      end

      if nargout > 2, relDiffs( iter ) = relDiff; end

      if numel( tol ) > 0  &&  tol > 0
        if relDiff < tol, break; end
      end
    end

    if restart == true && innerProd( Dgy, x - lastX ) > 0
      % Restart (kill momentum) when trajectory and -gradient form oblique angles
      k = 0;
      v = x;
    else
      k = k + 1;
      v = x + (1/theta) * ( x - lastX );
    end

  end

  if nargout > 1, objectiveValues = objectiveValues( 1 : iter ); end
  if nargout > 2, relDiffs = relDiffs( 1 : iter ); end

  xStar = x;
end

