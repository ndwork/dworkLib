
function [ xStar, oValues, relDiffs ] = stochasticProxGrad( x0, stepSize, gGrad, proxth, varargin )
  % Implements the stochastic proximal gradient method as detailed in
  % "Convergence of stochastic proximal gradient algorithm" by Rosasco et al.
  % for a function where the gradient has the form of a summation:  gGrad(x) = sum_{i=1}^N gGrad_i(x)
  %
  % [ out, oValues, relDiffs ] = stochasticProxGrad( x0, gGrad, proxth [, 'g', g, 'h', h, ...
  %   'nEpochs', nEpochs, 'nStoch', nStoch, 'verbose', true/false ] )
  %
  % Inputs:
  % x0 - initial guess for iterations
  % stepSize - scalar specifying the step size of the gradient descent step
  % gGrad - a function handle for a function of the prototype gGrad( x, indxs )
  %   where x is the current value of the input and indxs is a set of indices
  %   to use during the calculation
  %   If indxs is empty, then gGrad is to return the number of elements in the summation.
  % proxth - the proximal operator of the h function (with parameter t);
  %     two inputs: the vector and the scalar value of the parameter t
  %
  % Optional Inputs:
  % g - a handle to the g function.  This is needed to calculate the objective values.
  % h - a handle to the h function.  This is needed to calculate the objective values.
  % nEpochs - a scalar specifying the number of epochs to run (default is 10)
  % nStoch - the number of elements to use to estimate the gradient
  % verbose - true/false (default) to display informative statements
  %
  % Written by Nicholas Dwork - Copyright 2021
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addParameter( 'g', [] );
  p.addParameter( 'h', [] );
  p.addParameter( 'nEpochs', 10, @ispositive );
  p.addParameter( 'nStoch', 10, @ispositive );
  p.addParameter( 'saveEvery', 10, @ispositive );
  p.addParameter( 'verbose', false, @islogical );
  p.parse( varargin{:} );
  g = p.Results.g;
  h = p.Results.g;
  nEpochs = p.Results.nEpochs;
  nStoch = p.Results.nStoch;
  saveEvery = p.Results.saveEvery;
  verbose = p.Results.verbose;

  if nargin < 1
    disp( 'Usage:  [ out, oValues, relDiffs ] = stochasticProxGrad( x0, gGrad, proxth [, ' );
    disp( '  ''g'', g, ''h'', h, ''nEpochs'', nEpochs, ''nStoch'', nStoch, ' );
    disp( '  ''verbose'', true/false ] )' );
    if nargout > 0, xStar = []; end
    if nargout > 1, oValues = []; end
    if nargout > 2, relDiffs = []; end
    return
  end

  nGradTerms = gGrad( x0, [] );

  nSegs = ceil( nGradTerms / nStoch );  

  if nargout > 1, oValues = zeros( nEpochs * nSegs, 1 ); end
  if nargout > 2, relDiffs = zeros( nEpochs * nSegs, 1 ); end

  x = x0;  clear x0;
  iter = 0;
  for epochIndx = 1 : nEpochs
    indxs = randperm( nGradTerms );

    for segIndx = 1 : nSegs
      iter = iter + 1;
      subIndxs = indxs( (segIndx - 1 ) * nStoch + 1 : min( segIndx * nStoch, nGradTerms ) );

      Dgx = gGrad( x, subIndxs ) * ( nGradTerms / numel( subIndxs ) ) ;
      y = x - stepSize * Dgx;

      if nargout > 2, lastX = x; end
      x = proxth( y, stepSize );

      verboseStr = [ 'Epoch: ', indx2str( epochIndx, nEpochs ), ' of ', num2str(nEpochs), ',  ', ...
                     'Seg: ', indx2str( segIndx, nSegs ), ' of ', num2str(nSegs) ];

      if mod( iter-1, saveEvery ) == 0
        if nargout > 1
          objValue = g( x ) + h( x );
          oValues( ( epochIndx - 1 ) * nGradTerms + segIndx ) = objValue;
          verboseStr = [ verboseStr, ',  Obj Value: ', num2str( objValue ) ];   %#ok<AGROW>
        end

        if nargout > 2
          relDiff = norm( x(:) - lastX(:) ) / norm( lastX(:) );
          relDiffs( ( epochIndx - 1 ) * nGradTerms + segIndx ) = relDiff;
          verboseStr = [ verboseStr, ',  Rel Diff: ', num2str( relDiff ) ];   %#ok<AGROW>
        end
      end

      if verbose == true, disp( verboseStr ); end
    end

  end

  xStar = x;
end

