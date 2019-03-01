
function xStar = projSubgrad( x, gGrad, proj, varargin )
  % xStar = projSubgrad( x, gGrad, proj [, 'N', N, 't', t ] )
  %
  % This function implements the projected subgradient method
  %
  % Inputs:
  % x - the starting point
  % gGrad - a function handle representing the gradient function of g;
  %     input: the point to evaluation, output: the gradient vector
  % proj - a funciton handle representing the projection function
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
  p.addParameter( 'N', 100, @(x) isnumeric(x) && x>0 );
  p.addParameter( 't', 1, @(x) isnumeric(x) && x>0 );
  p.parse( varargin{:} );
  N = p.Results.N;
  t = p.Results.t;

  for n = 1:N
    x = x - t * gGrad( x );  % Gradient update
    x = proj( x );  % projection
  end

  xStar = x;
end
