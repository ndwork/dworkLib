
function out = proxShiftScale( f, x, a, b )
  % Computes the proximal operator of f( a x - b )
  %
  % Inputs:
  % f - a function handle to the proximal operator
  % x - the domain value for proximal operator evaluation
  % a - (scalar) the scaling of the domain variable
  % b - the shifting of the domain variable
  %     (scalar or array the same size as x)
  %
  % Written by Nicholas Dwork, Copyright 2019
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  out = 1/a * ( f( a*x + b, a*a ) - b );
end
