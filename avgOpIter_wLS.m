
function [x,objValues,alphasUsed] = avgOpIter_wLS( x0, S, varargin )
  % Implements an averaged opterator iteration with line search.
  % See "Line Search for Averaged Operator Iteration" by Gisellson et al. (2016)
  %
  % x = avgOpIter( x0, S [, 'alpha_bar', alpha_bar, 'maxIter', maxIter ] )
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
  p.addParameter( 'alpha_bar', 0.5, @(x) x>0 && x<1 );
  p.addParameter( 'doLineSearch', true, @islogical );
  %p.addParameter( 'doLineSearchTest', [], @islogical );
  p.addParameter( 'doLineSearchTest', false, @islogical );
  p.addParameter( 'N', 100, @ispositive );
  p.addParameter( 'objFunction', [] );
  p.addParameter( 'printEvery', 1, @ispositive );
  p.addParameter( 'verbose', false );
  p.parse( varargin{:} );
  alpha_bar = p.Results.alpha_bar;
  doLineSearch = p.Results.doLineSearch;
  doLineSearchTest = p.Results.doLineSearchTest;
  N = p.Results.N;
  objFunction = p.Results.objFunction;
  printEvery = p.Results.printEvery;
  verbose = p.Results.verbose;

  if numel( doLineSearchTest ) == 0, doLineSearchTest = doLineSearch; end
  if doLineSearchTest == true && doLineSearch == false
    error( 'Cannot run line search test without line search' );
  end

  if nargout > 1
    if numel( objFunction ) == 0
      error( 'Must specify an objective function to return objective values' );
    end
    objValues = zeros( N, 1 );
  end

  if nargout > 2
    alphasUsed = zeros( N, 1 );
  end
  
  %%% parameters
  eps = 0.01; % eps for (1 - eps) || rbar_k || in linesearch
  epsHat = 0.05;  % eps for line search test
  alpha0 = 2;
  %alpha_change = 1/1.4; % factor for change in alpha during linesearch
  alpha_change = 0.9; % factor for change in alpha during linesearch

  %alpha0 * ( alpha_change^k ) < alpha_bar
  %alpha0 / alpha_bar < 1 / ( alpha_change^k )
  %log( alpha0 / alpha_bar ) < -k log( alpha_change )
  %k > -log( alpha0 / alpha_bar ) / log( alpha_change )

  k = ceil( -log( alpha0 / alpha_bar ) / log( alpha_change ) );
  alphas = alpha0 .* ( alpha_change.^(0:k) );
  alphas(end) = alpha_bar;

  if isnumeric( S )
    S_fun = @(x) S * x;
  else
    S_fun = @(x) S(x);
  end

  x = x0;
  rk = S_fun(x) - x;
  normRk = sqrt( real( dotP( rk, rk ) ) );

  nAlphas = numel( alphas );
  normRks = zeros( nAlphas, 1 );
  xs = cell( 1, nAlphas );
  rks = cell( 1, nAlphas );

  for optIter = 1 : N
    %lastRk = rk;

    if doLineSearch == true

      parfor alphaIndx = 1 : nAlphas
        alpha = alphas( alphaIndx );
        xAlpha = x + alpha * rk;
        xs{alphaIndx} = xAlpha;
        rkAlpha = S_fun( xAlpha ) - xAlpha;   %#ok<PFBNS>
        rks{alphaIndx} = rkAlpha;
        normRks( alphaIndx ) = sqrt( real( dotP( rkAlpha, rkAlpha ) ) );
      end

      normRk_bar = normRks(end);
      normRks(end) = 0;
      bestAlphaIndx = find( normRks <= ( 1-eps ) * normRk_bar, 1 );

      alphaUsed = alphas( bestAlphaIndx );
      x = xs{ bestAlphaIndx };
      rk = rks{ bestAlphaIndx };

    else

      alphaUsed = alpha_bar;
      x = x + alpha_bar * rk;
      rk = S_fun( x ) - x;
    end

    if nargout > 1 || ( numel(objFunction) > 0 && verbose == true )
      objValue = objFunction( x );
    end

    if nargout > 1
      objValues( optIter ) = objValue;
    end
    if nargout > 2
      alphasUsed( optIter ) = alphaUsed;
    end

    if verbose == true && mod( optIter, printEvery ) == 0
      outStr = [ 'avgOpIter: Completed ', indx2str(optIter,N), ' of ', num2str(N), ...
        '   alpha: ', num2str( alphaUsed ) ];
      if numel( objFunction ) > 0
        outStr = [ outStr, '  objective: ', num2str( objValue ) ];   %#ok<AGROW>
      end
      disp( outStr );
    end

    if doLineSearchTest == true
      doLineSearch = false;
      lastNormRk = normRk;
      normRk = sqrt( real( dotP( rk, rk ) ) );
      if normRk * lastNormRk == 0, continue; end
  
      %if ( real( dotP( rk, lastRk ) ) / ( normRk * lastNormRk ) ) > ( 1 - epsHat )
      if ( alphaUsed ~= alpha_bar )  ||  ( normRk / lastNormRk < 1 - epsHat )  % My own heuristic
        % The values of S are changing quickly in this region
        doLineSearch = true;
      end
    end

  end

end
