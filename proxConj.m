
function out = proxConj( proxf, in, sigma )
  % computes the proximal operator of the scaled conjugate function of f
  %
  % Inputs:
  % in - the vector to compute the proximal operator of
  % proxf - a function handle to the proximal operator of f
  % sigma - the scaling of the conjugate function
  %
  % Written by Nicholas Dwork, Copyright 2024
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  out = in - sigma * proxf( in / sigma, 1 / sigma );
end
