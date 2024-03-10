
function [x,objValues] = avgOpIter( x0, S, varargin )
  % Implements and averaged opterator iteration.  See "Line Search for Averaged
  % Operator Iteration" by Gisellson et al. (2016)
  %
  % x = avgOpIter( x0, S [, 'alpha', alpha, 'N', N ] )
  %
  % Inputs:
  % x0 - the initial guess
  % S - either a matrix or a a function handle that is the non-expansive operator
  %
  % Optional Inputs:
  % alpha - the scalar for the combination of the averaged operator iteration (default is 0.5)
  % N - the number of iterations to run (default is 100)
  % objFunction - a function handle to the objective function
  % printEvery - print verbose statements every printEvery iterations
  %
  % Outputs:
  % x - the value of the domain variable after all iterations are complete
  % objValues - the value of objFunction evaluated on x at each iteration
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
  p.addParameter( 'printEvery', 1, @ispositive );
  p.addParameter( 'verbose', false );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  N = p.Results.N;
  objFunction = p.Results.objFunction;
  printEvery = p.Results.printEvery;
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

    if verbose == true && mod( optIter, printEvery ) == 0
      outStr = [ 'avgOpIter: Completed ', indx2str(optIter,N), ' of ', num2str(N) ];
      if numel( objFunction ) > 0
        outStr = [ outStr, '  objective: ', num2str( objValue ) ];   %#ok<AGROW>
      end
      disp( outStr );
    end

    objValues(optIter) = objFunction( x );
  end

end
