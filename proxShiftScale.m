
function out = proxShiftScale( proxf, x, a, b )
  % Computes the proximal operator of f( a x + b )
  %
  % Inputs:
  % proxf - a function handle to the proximal operator of f
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

  out = (1/a) * ( proxf( a*x + b, a*a ) - b );
end
