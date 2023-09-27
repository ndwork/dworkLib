
function [xStar,objValues,relDiffs] = pdhg( x, proxf, proxgConj, tau, varargin )
  % [xStar,objValues,relDiffs] = pdhg( x, proxf, proxgConj, tau [, ...
  %   'A', A, 'f', f, 'g', g, 'N', N, 'normA', normA, 'sigma', sigma, ...
  %    'lambda', lambda, 'printEvery', printEvery, 'theta', theta, ...
  %    'tol', tol, 'verbose', verbose, 'z', z ] )
  %
  % Implements the Primal Dual Hybrid Gradient (Chambolle-Pock) method that
  % solves problems of the form:  minimize f( x ) + g( A x )
  %
  % Inputs:
  % x - initial guess
  % proxf - a function handle for the proximal operator of
  % proxgConj - a function handle for the proximal operator of the conjugate of g
  %
  % Optional Inputs:
  % A - if A is not provided, it is assumed to be the identity
  % f - to determine the objective values, f must be provided
  % g - to determine the objective values, g must be provided
  % N - the number of iterations that ADMM will perform (default is 100)
  % normA - the matrix induced 2-norm of A.  Could be determined with norm or, for
  %   large sparse matries, estimated with normest or powerIteration.
  % lambda - relaxation parameter
  % theta - acceleration parameter
  % verbose - true or false
  % z - initial value of dual variable
  %
  % Outputs:
  % xStar - the optimal point
  %
  % Optional Outputs:
  % objValues - a 1D array containing the objective value of each iteration
  % relDiffs - a 1D array containing the relative difference of each iteration
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp( 'Usage;  ' );
    disp( '  [xStar,objValues] = pdhg( x, proxf, proxgConj, tau [, ... ');
    disp( '  ''A'', A, ''f'', f, ''g'', g, ''N'', N, ''normA'', normA, ''sigma'', sigma, ...' );
    disp( '  ''lambda'', lambda, ''printEvery'', printEvery, ''theta'', theta, ...' );
    disp( '  ''tol'', tol, ''verbose'', verbose, ''z'', z ] )' );
    if nargout > 0, xStar = []; end
    if nargout > 1, objValues = []; end
    return;
  end

  defaultTol = 1d-4;

  p = inputParser;
  p.addParameter( 'A', [] );
  p.addParameter( 'f', [] );
  p.addParameter( 'g', [] );
  p.addParameter( 'N', 100, @ispositive );
  p.addParameter( 'normA', [], @(x) x>0 );
  p.addParameter( 'sigma', [], @ispositive );
  p.addParameter( 'theta', 1, @(x) x >= 0 && x <= 1 );
  p.addParameter( 'tol', defaultTol, @(x) numel(x) == 0 || ispositive(x) );
  p.addParameter( 'lambda', 1, @(x) x >= 0 && x <= 2 );
  p.addParameter( 'verbose', false, @(x) islogical(x) || x==1 || x==0 );
  p.addParameter( 'printEvery', 1, @ispositive );
  p.addParameter( 'z', [], @isnumeric );
  p.parse( varargin{:} );
  A = p.Results.A;
  f = p.Results.f;
  g = p.Results.g;
  N = p.Results.N;
  normA = p.Results.normA;
  sigma = p.Results.sigma;
  theta = p.Results.theta;
  tol = p.Results.tol;
  lambda = p.Results.lambda;
  printEvery = p.Results.printEvery;
  verbose = p.Results.verbose;
  z = p.Results.z;
  
  if numel( A ) == 0
    applyA = @(x) x;
    applyAT = @(x) x;
  elseif isnumeric( A )
    applyA = @(x) A * x;
    applyAT = @(y) A' * y;
  else
    applyA = @(x) A( x, 'notransp' );
    applyAT = @(x) A( x, 'transp' );
  end

  if nargout > 1,  objValues = zeros( N, 1 ); end
  if nargout > 2,  relDiffs = zeros( N, 1 ); end

  if numel( sigma ) == 0  && numel( A ) > 0
    if numel( normA ) == 0
      error( 'If an A is supplied, you must supply sigma or normA' );
    end
    if normA == 0, return; end
    sigma = ( 0.99 / normA^2 ) / tau;
  end

  calculateRelDiffs = false;
  if nargout > 2  ||  ( numel( tol ) > 0  &&  tol < Inf )
    calculateRelDiffs = true;
  end

  maxAbsX = max( abs( x(:) ) );
  xBar = zeros( size (x) );
  if numel( z ) == 0, z = 0; end
  
  optIter = 0;
  relDiff = Inf;
  while optIter < N
    optIter = optIter + 1;

    zTmp = z + sigma * applyA( xBar );
    z = proxgConj( zTmp, sigma );

    lastX = x;
    lastMaxAbsX = maxAbsX;

    ATz = applyAT( z );
    xTmp = x - tau * ATz;
    x = proxf( xTmp, tau );

    xBar = x + theta * ( x - lastX );

    maxAbsX = max( abs( x(:) ) );

    if nargout > 1
      objValues( optIter ) = f( x ) + g( applyA( x ) );
    end
    if calculateRelDiffs == true && optIter > 1
      normLastX = norm( lastX(:) );
      normDiffX = norm( x(:) - lastX(:) );
      if normLastX == 0
        if normDiffX == 0
          relDiff = 0;
        else
          relDiff = normDiffX / norm( x(:) );
        end
      else
        relDiff = normDiffX / normLastX;
      end
    end
    if nargout > 2, relDiffs( optIter ) = relDiff; end

    if verbose == true
      verboseStr = [ 'pdhg: iteration ', indx2str(optIter,N), ' of ', num2str(N) ];
      if nargout > 1
        verboseStr = [ verboseStr, ',  objective value: ', num2str( objValues( optIter ) ) ];   %#ok<AGROW>
      end
      if calculateRelDiffs == true
        verboseStr = [ verboseStr, ',  relative diff: ', num2str( relDiff ) ];   %#ok<AGROW>
      end
      if mod( optIter, printEvery ) == 0  ||  optIter == 1, disp( verboseStr ); end
    end

    if optIter > 1  &&  numel( tol ) > 0  &&  tol > 0  &&  tol < Inf
      if relDiff < tol, break; end
    end
    if maxAbsX == 0  && lastMaxAbsX == 0 && optIter > 1
      break;
    end
  end

  if nargout > 1, objValues = objValues( 1 : optIter ); end
  if nargout > 2, relDiffs = relDiffs( 1 : optIter ); end

  xStar = x;
end
