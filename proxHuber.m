
function out = proxHuber( x, t )
  % out = proxHuber( x, t, mu )
  %
  % Calculates the proximal operator of f(x) = t * Huber( x, mu ).
  % If x is an array, proxHuber operates on each element of the array
  %   individually.
  %
  % Inputs:
  % x - an
  % t - a scalar
  % mu - the Huber penalty parameter
  %
  % Written by Nicholas Dwork - Copyright 2021
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  if nargin < 1
    disp( 'Usage:  out = proxHuber( x, t, mu )' );
    out = [];
    return
  end

  out = ( t .* softThresh( x, t + mu ) + mu * x ) ./( t + mu );

end
