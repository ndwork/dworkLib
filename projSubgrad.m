
function [xStar,objValues] = projSubgrad( x, gGrad, proj, varargin )
  % [xStar,objValues] = projSubgrad( x, gGrad, proj [, 'g', g, 'N', N, 't', t ] )
  %
  % This function implements the projected subgradient method
  %
  % Inputs:
  % x - the starting point
  % gGrad - a function handle that returns a (sub)gradient of g
  %   input: the point to evaluation, output: the (sub)gradient vector
  % proj - a function handle representing the projection function
  %
  % Optional Inputs:
  % alpha - line search parameter for Armijo line search
  % beta - line search parameter for step size reduction (default is 0.8)
  % N - the number of iterations that will be performed (default is 100)
  % t - step size (default is 1)
  %
  % Outputs:
  % xStar - the optimal point
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addParameter( 'alpha', 0.3, @(x) x > 0 && x < 0.5 );
  p.addParameter( 'beta', 0.8, @(x) x > 0 && x < 1 );
  p.addParameter( 'g', [] );
  p.addParameter( 'N', 100, @(x) isnumeric(x) && x>0 );
  p.addParameter( 'nMaxLineSearchIter', 100, @ispositive );
  p.addParameter( 'printEvery', 1, @isnumeric );
  p.addParameter( 't', 1, @(x) isnumeric(x) && x>0 );
  p.addParameter( 'tau', 1.1, @(x) x >= 1 );
  p.addParameter( 'tThresh', 1d-18 );
  p.addParameter( 'useLineSearch', true, @islogical );
  p.addParameter( 'verbose', 0, @(x) isnumeric(x) || islogical(x) );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  beta = p.Results.beta;
  g = p.Results.g;
  N = p.Results.N;
  nMaxLineSearchIter = p.Results.nMaxLineSearchIter;
  printEvery = p.Results.printEvery;
  t = p.Results.t;
  tau = p.Results.tau;
  tThresh = p.Results.tThresh;
  useLineSearch = p.Results.useLineSearch;
  verbose = p.Results.verbose;

  if nargin < 1
    xStar = [];  objValues = [];
    disp( 'Usage:  [xStar,objValues] = projSubgrad( x, gGrad, proj [, ''g'', g,' );
    disp( '          ''N'', N, ''t'', t ] )' );
    return
  end

  if useLineSearch == true  &&  numel( g ) == 0
    error( 'Must supply g to use line search' );
  end

  if nargout > 1  ||  useLineSearch == true
    gx = g( x );
  end

  if nargout > 1
    objValues = zeros( N+1, 1 );
    objValues( 1 ) = gx;
  end

  for n = 1:N
    if verbose ~= 0  &&  mod( n, printEvery ) == 0
      disp([ 'Working on iteration ', num2str(n), ' of ', num2str(N) ]);
    end

    gGradX = gGrad( x );
    normGradX = norm( gGradX(:) );

    if useLineSearch == true
      lineBreakThresh = gx - alpha * t * normGradX^2;
      for lineSearchIter = 1 : nMaxLineSearchIter
        xNew = x - t * gGradX;
        gxNew = g( xNew );
        if gxNew < lineBreakThresh || t < tThresh
          break;
        end
        t = beta * t;
        if t < 1d-18
          break;
        end
      end
      t = tau * t;
      x = xNew;
      gx = gxNew;
    else
      if nargout > 1, gx = g( x ); end
      x = x - t * gGradX;
    end

    x = proj( x );

    if nargout > 1, objValues( n+1 ) = gx; end
  end

  xStar = x;
end
