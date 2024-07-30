
function out = LpNorms( A, p, dim )
  % out = LpNorms( A [, p, dim ] )
  %
  % Calculates the Lp norm of A along dimension dim
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
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp( 'Usage:  out = LpNorms( A [, p, dim ] )' );
    return
  end
  
  if nargin < 3, dim = ndims(A); end
  if nargin < 2, p = 2; end

  absA = abs( A );

  if mod( p, 1 ) == 0
    % p is an integer
    tmp = absA;
    for i=2:p
      tmp = tmp .* absA;
    end    
  else
    % p is not an integer
    tmp = absA.^p;
  end

  tmp = sum( tmp, dim );
  out = tmp .^ (1/p);
end
