
function out = numericalDiff( g, x, varargin )
  % out = numericalDiff( g, x [, 'dx', dx ] );
  %
  % Numerically estimates the gradient of function g
  %
  % Inputs:
  % g - function handle that accepts a variable with number of elements of x
  % x - the input where the derivative is to be estimated
  %
  % Outputs:
  % the estimate of the gradient
  %
  % Written by Nicholas Dwork - Copyright 2023
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  if nargin < 2
    disp( 'Usage:  out = numericalDiff( g, x [, ''dx'', dx ] );')
    out = [];
    return;
  end

  p = inputParser;
  p.addParameter( 'dx', 0.01, @ispositive );
  p.parse( varargin{:} );
  dx = p.Results.dx;

  out = zeros( size( x ) );
  xH = x;
  xL = x;
  for i = 1 : numel( x )
    xH(i) = x(i) + dx;
    gxH = g( xH );
    xL(i) = x(i) - dx;
    gxL = g( xL );
    out(i) = ( gxH - gxL ) / ( 2 * dx );

    % reset xH and xL
    xH(i) = x(i);  xL(i) = x(i);
  end

end
