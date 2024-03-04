function x = avgOpIter( x0, S, varargin )
  % Implements and averaged opterator iteration.  See "Line Search for Averaged
  % Operator Iteration" by Gisellson et al. (2016)
  %
  % x = avgOpIter( x0, S [, 'alpha', alpha, 'maxIter', maxIter ] )
  %
  % Inputs:
  % x0 - the initial guess
  % S - either a matrix or a a function handle that is the non-expansive operator
  %
  % Written by Nicholas - Copyright 2024
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addParameter( 'alpha', 0.5, @(x) x>0 && x<1 );
  p.addParameter( 'maxIter', 100, @ispositive );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  maxIter = p.Results.maxIter;

  if ismatrix( S )
    x = x0;
    for i = 1 : maxIter
      x = ( 1 - alpha ) * x + alpha * S * x;
    end

  else
    % S is a function handle
  
    x = x0;
    for i = 1 : maxIter
      x = ( 1 - alpha ) * x + alpha * S( x );
    end
  end

end
