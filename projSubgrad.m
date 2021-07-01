
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
  p.addParameter( 'g', [] );
  p.addParameter( 'N', 100, @(x) isnumeric(x) && x>0 );
  p.addParameter( 't', 1, @(x) isnumeric(x) && x>0 );
  p.addParameter( 'verbose', 0, @(x) isnumeric(x) || islogical(x) );
  p.parse( varargin{:} );
  g = p.Results.g;
  N = p.Results.N;
  t = p.Results.t;
  verbose = p.Results.verbose;

  if nargin < 1
    xStar = [];  objValues = [];
    disp( 'Usage:  [xStar,objValues] = projSubgrad( x, gGrad, proj [, ''g'', g,' );
    disp( '          ''N'', N, ''t'', t ] )' );
    return
  end

  if nargout > 1
    objValues = zeros( N+1, 1 );
    objValues( 1 ) = g( x );
  end

  for n = 1:N
    if verbose ~= 0
      disp([ 'Working on iteration ', num2str(n), ' of ', num2str(N) ]);
    end

    x = x - t * gGrad( x );  % Gradient update
    x = proj( x );  % projection

    if nargout > 1, objValues( n+1 ) = g( x ); end
  end

  xStar = x;
end
