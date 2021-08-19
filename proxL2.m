
function out = proxL2( x, t )
  % out = proxL2( x, t )
  %
  % Calculates the proximal operator of f(x) = t * L2( x )
  %
  % Inputs:
  % x - a 1D array
  % t - a scalar
  %
  % Written by Nicholas Dwork - Copyright 2020
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  nx = norm( x(:) );

  if nx <= t
    out = 0;
    return;
  end

  out = ( 1 - t / nx ) * x;
end
