
function [x,objValues] = avgOpIter( x0, S, varargin )
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
  p.addParameter( 'N', 100, @ispositive );
  p.addParameter( 'objFunction', [] );
  p.addParameter( 'verbose', false );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  N = p.Results.N;
  objFunction = p.Results.objFunction;
  verbose = p.Results.verbose;

  if nargout > 1
    if numel( objFunction ) == 0
      error( 'Must specify an objective function to return objective values' );
    end
    objValues = zeros( N, 1 );
  end

  x = x0;
  for optIter = 1 : N
    if isa( S, 'function_handle' )
      x = ( 1 - alpha ) * x + alpha * S( x );
    else
      x = ( 1 - alpha ) * x + alpha * S * x;
    end

    if nargout > 1 || ( numel(objFunction) > 0 && verbose == true )
      objValue = objFunction( x );
    end

    if nargout > 1
      objValues( optIter ) = objValue;
    end

    if verbose == true
      outStr = [ 'avgOpIter: Completed ', indx2str(optIter,N), ' of ', num2str(N) ];
      if numel( objFunction ) > 0
        outStr = [ outStr, '  objective: ', num2str( objValue ) ];   %#ok<AGROW>
      end
      disp( outStr );
    end

    objValues(optIter) = objFunction( x );
  end

end
