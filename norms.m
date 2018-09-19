
function out = norms( A, varargin )
  % out = norms( A [, p, dim ] )
  %
  % Calculates the p norm of A along dimension dim
  %
  % Inputs:
  % A - a matrix of dimension <= dim
  %
  % Optional Inputs:
  % p - Compute the Lp norm (default is 2)
  % dim - (default is last dimension)
  %
  % Written by Nicholas Dwork - Copyright 2018
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  par = inputParser;
  par.addOptional( 'p', 2, @isnumeric );
  par.addOptional( 'dim', ndims(A), @isnumeric );
  par.parse( varargin{:} );
  p = par.Results.p;
  dim = par.Results.dim;

  if mod( p, 1 ) == 0
    % p is an integer
    tmp = abs(A);
    for i=2:p
      tmp = tmp .* abs(A);
    end    
  else
    % p is not an integer
    tmp = abs(A).^p;
  end

  tmp = sum( tmp, dim );
  out = tmp .^ (1/p);
end
