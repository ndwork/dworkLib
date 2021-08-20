
function out = proxHuber( x, t, mu, b )
  % out = proxHuber( x, t, mu, b )
  %
  % Calculates the proximal operator of f(x) = t * Huber( x - b, mu ).
  % If x is an array, proxHuber operates on each element of the array
  %   individually.
  %
  % Inputs:
  % x - an
  % t - a scalar
  % mu - the Huber penalty parameter
  % b - a translation parameter
  %
  % Written by Nicholas Dwork - Copyright 2021
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  if nargin < 1
    disp( 'Usage:  out = proxHuber( x [, t, mu, b ] )' );
    out = [];
    return
  end

  if nargin < 2, t = 1; end
  if nargin < 3, mu = 1; end

  if nargin < 4  % b = 0;
    out = ( t .* softThresh( x, t + mu ) + mu * x ) ./( t + mu );

  else
    out = ( t .* softThresh( x - b, t + mu ) + mu * x ) ./ ( t + mu ) + b;

  end

end
