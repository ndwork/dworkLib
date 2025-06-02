
function out = projectOntoBall( x, r, p )
  % Project onto a ball of radius r with respect to the Lp norm.
  %
  % out = projectOntoBall( x, r [, p ] )
  %
  % Inputs -
  % x - 1D array that is the vector to project
  % r - the radius of the ball to project onto
  %
  % Optional inputs:
  % p - either a positive real scalar, Inf, -Inf, or 'fro'
  %
  % Outputs:
  % out - a 1D array that is the projected vector
  %
  % Written by Nicholas Dwork, Copyright 2024
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 3, p = 2; end

  xNorm = norm( x, p );

  if xNorm <= r
    out = x;
  else
    out = x * ( r / xNorm );
  end
end
