
function out = projectOntoProbSimplex( y )
  % out = projectOntoProbSimplex( y )
  %
  % This function performs a Euclidean projection of y onto the probability simplex.
  % That is, out is the nearest point to y such that out is non-negative and it
  % sums to 1.
  %
  % This function is based on "Projection onto the probability simplex: An 
  % efficient algorithm with a simple proof, and an application" by Wang and
  % Carreira-Perpinan.
  %
  % Written by Nicholas Dwork, Copyright 2021
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  if nargin < 1
    disp( 'Usage:  out = projectOntoProbSimplex( y )' );
    out = [];  return;
  end

  u = sort( y, 'descend' );
  j = ( 1 : numel(y) )';

  rho = find( ( u + ( 1 ./ j ) .* ( 1 - cumsum( u(:) ) ) ) > 0, 1, 'last' );

  lambda = (1/rho) * ( 1 - sum( u(1:rho) ) );

  out = max( y + lambda, 0 );
end
