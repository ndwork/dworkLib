
function [xStar,objValues,metricValues] = pdhgWLS( x, proxf, proxgConj, varargin )
  % [xStar,objValues,metricValues] = pdhgWLS( x, proxf, proxgConj [, ...
  %   'dsc', dsc, 'dscThresh', dscThresh, 'N', N, 'A', A, 'beta', beta, ...
  %   'f', f, 'g', g, 'mu', mu, 'P', P, 'tau', tau, 'theta', theta, 'y', y, 'verbose', verbose ] )
  %
  % Implements Primal-Dual Hybrid graident method (Chambolle-Pock) with line search
  % based on "A First-Order Primal-Dual Algorithm with Linesearch" by Malitsky and Pock
  %
  % minimizes f( x ) + g( A x )
  %
  % OR, with preconditioning,
  % minimizes f( y ) + g( A P^{-1} y ).  Then y = P x  <==>  x = P^{-1} y;
  %
  % Inputs:
  % x - initial guess
  %
  % Optional Inputs:
  % A - if A is not provided, it is assumed to be the identity
  % beta - line search parameter
  % dsc - use dynamic stopping criteria (true / false) or a function handle that is the dynamic
  %       stopping criteria function.  This function must accept x and lastX as inputs and output
  %       a scalar value.
  % dscThresh - when dsc < dscThresh, iterations stop
  % f - to determine the objective values
  %     If if is empty, it is assumed it is the 0 function
  % g - to determine the objective values, g must be provided
  % metrics - either a function handle or a cell array of function handles for functions that
  %           should be run on every iteration and the values reported when verbose is true
  % N - the maximum number of iterations that PDHG will perform (default is 1000)
  % tau - the step size parameter that gets altered with line search (default is 1)
  % y - the initial values of y in the PDHG iterations
  %
  % Outputs:
  % xStar - the optimal point
  %
  % Optional Outputs:
  % objValues - a 1D array containing the objective value of each iteration
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 3
    disp( 'Usage:   [xStar,objValues,metricValues] = pdhgWLS( x, proxf, proxgConj [, ... ' );
  	disp( '           ''dsc'', dsc, ''dscThresh'', dscThresh, ''N'', N, ''A'', A, ''beta'', beta, ... ' );
    disp( '           ''f'', f, ''g'', g, ''metrics'', metrics, ''mu'', mu, ''P'', P, ... ' );
    disp( '           ''tau'', tau, ''theta'', theta, ''y'', y, ''verbose'', verbose ] ) ' );
    if nargout > 0, xStar = []; end
    if nargout > 1, objValues = []; end
    return;
  end
  
  defaultN = 1000;

  p = inputParser;
  p.addParameter( 'A', [] );
  p.addParameter( 'beta', 1, @ispositive );
  p.addParameter( 'delta', 0.99, @(x) x>0 && x<1 );
  p.addParameter( 'doCheckAdjoint', false, @(x) islogical(x) || x == 1 || x == 0 );
  p.addParameter( 'dsc', false );
  p.addParameter( 'dscThresh', 1d-8, @ispositive );
  p.addParameter( 'f', [] );
  p.addParameter( 'g', [] );
  p.addParameter( 'innerProd', [] );
  p.addParameter( 'metricNames', [] );
  p.addParameter( 'metrics', [] );
  p.addParameter( 'mu', 0.8, @(x) x>0 && x<1 );
  p.addParameter( 'N', defaultN, @ispositive );
  p.addParameter( 'P', [] );
  p.addParameter( 'printEvery', 1, @ispositive );
  p.addParameter( 'saveDir', './', @(x) true );
  p.addParameter( 'saveEvery', [], @ispositive );
  p.addParameter( 'saveFun', [] );
  p.addParameter( 'tau', [], @ispositive );
  p.addParameter( 'theta', 1, @ispositive );
  p.addParameter( 'y', [], @isnumeric );
  p.addParameter( 'verbose', false, @(x) islogical(x) || x == 1 || x == 0 );
  p.parse( varargin{:} );
  A = p.Results.A;
  beta = p.Results.beta;
  delta = p.Results.delta;
  doCheckAdjoint = p.Results.doCheckAdjoint;
  dsc = p.Results.dsc;
  dscThresh = p.Results.dscThresh;
  f = p.Results.f;
  g = p.Results.g;
  innerProd = p.Results.innerProd;
  metricNames = p.Results.metricNames;
  metrics = p.Results.metrics;
  mu = p.Results.mu;
  N = p.Results.N;
  P = p.Results.P;
  printEvery = p.Results.printEvery;
  saveDir = p.Results.saveDir;
  saveEvery = p.Results.saveEvery;
  saveFun = p.Results.saveFun;
  tau = p.Results.tau;
  theta = p.Results.theta;
  y = p.Results.y;
  verbose = p.Results.verbose;

  if nargout > 1
    if numel( g ) == 0, error( 'Must supply g to calculate the objective values' ); end
  end

  if numel( tau ) == 0, tau = 1; end

  if numel( innerProd ) == 0
    innerProd = @(x,y) real( dotP( x, y ) );
  end

  if numel( N ) == 0, N = defaultN; end

  if numel( A ) == 0  &&  numel( P ) == 0
    applyA = @(x) x;
    applyAT = @(x) x;

  elseif numel( A ) > 0  &&  numel( P ) == 0
    if isnumeric( A )
      applyA = @(x) A * x;
      applyAT = @(y) A' * y;
    else
      applyA = @(x) A( x, 'notransp' );
      applyAT = @(x) A( x, 'transp' );
    end

  elseif numel( A ) == 0 && numel( P ) > 0
    if isnumeric( P )
      applyA = @(x) P \ x;
      applyAT = @(x) P' \ x;
    else
      applyA = @(x) P( x, 'inv' );
      applyAT = @(x) P( x, 'invTransp' );
    end

  elseif numel( A ) > 0 && numel( P ) > 0
    if ~isnumeric( A ) && ~isnumeric( P )
      applyA = @(x) A( P( x, 'inv' ) );
      applyAT = @(x) P( A( x, 'transp' ), 'invTransp' );
    elseif ~isnumeric( A ) && isnumeric( P )
      applyA = @(x) A( P \ x );
      applyAT = @(x) P' \ A( x, 'transp' );
    elseif isnumeric( A ) && ~isnumeric( P )
      applyA = @(x) A * P( x, 'inv' );
      applyAT = @(x) P( A' * x, 'invTransp' );    
    elseif isnumeric( A ) && isnumeric( P )
      applyA = @(x) A * ( P \ x );
      applyAT = @(x) P' \ ( A * x );
    end
  end

  if numel( y ) == 0
    y = applyA( x );
    y(:) = 0;
  end

  nMetrics = numel( metrics );
  if nMetrics == 1  &&  isa( metrics, "function_handle" )
    metrics = { metrics };
  end
  if isa( metricNames, "char" )
    metricNames = { metricNames };
  end
  if numel( metricNames ) > 0  &&  numel( metricNames ) ~= nMetrics
    error( 'Must either supply no metric names or the same number as metrics' );
  end

  relDiff = @(x,lastX) norm( x(:) - lastX(:) ) / norm( lastX(:) );
  if numel( dsc ) > 0  && ~isa( dsc, "function_handle" ) &&  dsc == true
    dsc = relDiff;
  end

  if ~exist( saveDir, 'dir' ), mkdir( saveDir ); end

  if doCheckAdjoint == true
    [adjointCheckPassed,adjCheckErr] = checkAdjoint( x, applyA, applyAT );
    if ~adjointCheckPassed, error([ 'checkAdjoint failed with error ', num2str(adjCheckErr) ]); end
  end

  if nargout > 1
    fx = 0;
    objValues = zeros( N+1, 1 );
  end
  if nargout > 2
    metricValues = zeros( N+1, nMetrics );
  else
    metricValues = zeros( nMetrics, 1 );
  end

  for optIter = 1 : N
    if nargout > 1
      if numel( f ) > 0, fx = f( x ); end
      gAx = g( applyA( x ) );
      objValues( optIter ) = fx + gAx;
    end

    if numel( saveEvery ) > 0  &&  ( mod( optIter, saveEvery ) == 0  ||  optIter == 1 )
      if numel( saveFun ) > 0
        saveFun( x, optIter, saveDir );
      else
        save( [ saveDir, '/pdhgWLS_save_', indx2str(optIter,N), '.mat' ], 'x' );
      end
    end

    if verbose == true
      dispStr = [];
      if nMetrics > 0
        for mIndx = 1 : nMetrics
          mValue = metrics{ mIndx }( x );
          if nargout > 2
            metricValues( optIter+1, mIndx ) = mValue;
          else
            metricValues( mIndx ) = mValue;
          end
        end
      end
      if mod( optIter, printEvery ) == 0 || optIter == 1
        dispStr = [ 'pdhgWLS: ', indx2str(optIter,N), ' of ', num2str(N) ];
        if nargout > 1
          dispStr = [ dispStr, ',  objective: ', num2str( objValues( optIter ),'%15.13f' ) ];   %#ok<AGROW>
        end
        if nMetrics > 0
          for mIndx = 1 : nMetrics
            if nargout > 2
              mValue = metricValues( optIter+1, mIndx );
            else
              mValue = metricValues( mIndx );
            end
            if numel( metricNames ) > 0
              dispStr = [ dispStr, ',  ', metricNames{mIndx}, ': ', num2str(mValue) ];   %#ok<AGROW>
            else
              dispStr = [ dispStr, ',  metric ', indx2str(mIndx,nMetrics), ': ', num2str(mValue) ];   %#ok<AGROW>
            end
          end
        end
      end
    end

    lastX = x;
    x = lastX - tau * applyAT( y );
    if numel( proxf ) > 0
      x = proxf( x, tau );
    end

    if isa( dsc, "function_handle" )
      dscValue = dsc( x, lastX );
      if verbose == true && numel( dispStr ) > 0
        dispStr = [ dispStr, ',  dscValue: ', num2str(dscValue) ];   %#ok<AGROW>
        disp( dispStr );
      end
      if dscValue < dscThresh  &&  optIter > 1, break; end
    else
      if verbose == true, disp( dispStr ); end
    end

    lastTau = tau;
    tau = tau * sqrt( 1 + theta );

    diffx = x - lastX;

    lastY = y;
    subIter = 0;
    while true
      subIter = subIter + 1;
      if verbose ~= false && printEvery == 1
        disp([ '     sub iteration: ', num2str( subIter ), '  tau: ', num2str(tau) ]);
      end

      theta = tau / lastTau;
      xBar = x + theta * diffx;

      betaTau = beta * tau;
      tmp = lastY + betaTau * applyA( xBar );
      y = proxgConj( tmp, betaTau );

      diffy = y - lastY;
      ATdiffy = applyAT( diffy );

      normATdiffy = sqrt( innerProd( ATdiffy(:), ATdiffy(:) ) );
      normDiffy = sqrt( innerProd( diffy(:), diffy(:) ) );

      if tau * sqrt( beta ) * normATdiffy <= delta * normDiffy
        break
      end
      tau = mu * tau;
    end
  end

  if numel( P ) > 0
    if isnumeric( P )
      x = P \ x;
    else
      x = P( x, 'inv' );
    end
  end

  xStar = x;
end
