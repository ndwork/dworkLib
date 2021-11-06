
function [ xStar, oValues, relDiffs ] = proxSVRG( x0, stepSize, gGrad, proxth, varargin )
  % Implements the accelerated stochastic proximal gradient method with accelerated
  % described in "Stochastic Proximal Gradient Descent with Acceleration Techniques"
  % by Nitanda
  % The gradient has the form of a summation:  gGrad(x) = sum_{i=1}^n gGrad_i(x)
  % Note that this differs from gGrad of the paper by a factor of n.
  %
  % [ out, oValues, relDiffs ] = proxSVRG( x0, gGrad, proxth [, 'g', g, 'h', h, ...
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
  % b - minibatch size
  % beta - acceleration factor.  This variable may have several forms:
  %   Empty: uses the default values (the same as FISTA)
  %   Scalar:  uses the same value for all iterations
  %   A function that accepts a positive integer k and returns a non-negative number
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
  p.addParameter( 'b', [], @ispositive );
  p.addParameter( 'beta', [] );
  p.addParameter( 'g', [] );
  p.addParameter( 'h', [] );
  p.addParameter( 'nS', 10, @ispositive );
  p.addParameter( 'nK', 10, @ispositive );
  p.addParameter( 'saveEvery', 10, @ispositive );
  p.addParameter( 'verbose', false, @islogical );
  p.parse( varargin{:} );
  b = p.Results.b;
  beta = p.Results.b;
  g = p.Results.g;
  h = p.Results.g;
  nS = p.Results.nS;
  nK = p.Results.nK;
  saveEvery = p.Results.saveEvery;
  verbose = p.Results.verbose;

  if nargin < 1
    disp( 'Usage:  [ out, oValues, relDiffs ] = proxSVRG( x0, gGrad, proxth [, ''b'', b, ' );
    disp( '  ''g'', g, ''h'', h, ''nS'', nS, ''nK'', nK, ''verbose'', true/false ] )' );
    if nargout > 0, xStar = []; end
    if nargout > 1, oValues = []; end
    if nargout > 2, relDiffs = []; end
    return
  end

  nGradTerms = gGrad( x0, [] ); 

  if nargout > 1, oValues = zeros( nS * nK, 1 ); end
  if nargout > 2, relDiffs = zeros( nS * nK, 1 ); end

  xTilde = x0;  clear x0;
  iter = 0;
  for s = 1 : nS

    vTilde = gGrad( x, 1:nGradTerms );
    x = xTilde;
    y = xTilde;

    for k = 1 : nK
      iter = iter + 1;

      indxs = randperm( nGradTerms );
      subIndxs = indxs( 1 : b );

      v = gGrad( y, subIndxs ) - gGrad( xTilde, subIndxs ) + vTilde;
      
      tmp = y - stepSize * v;
      lastX = x;
      x = proxth( tmp, stepSize );

      if numel( beta ) == 0
        thisBeta = k / ( k+3 );
      elseif isnumeric( beta )
        thisBeta = beta;
      else
        thisBeta = beta( k );
      end
      y = x + thisBeta * ( x - lastX );

      verboseStr = [ 's: ', indx2str( s, nS ), ' of ', num2str( nS ), ',  ', ...
                     'Seg: ', indx2str( k, nK ), ' of ', num2str( nK ) ];

      if mod( iter-1, saveEvery ) == 0
        if nargout > 1
          objValue = g( x ) + h( x );
          oValues( ( s - 1 ) * nGradTerms + k ) = objValue;
          verboseStr = [ verboseStr, ',  Obj Value: ', num2str( objValue ) ];   %#ok<AGROW>
        end

        if nargout > 2
          relDiff = norm( x(:) - lastX(:) ) / norm( lastX(:) );
          relDiffs( ( s - 1 ) * nGradTerms + k ) = relDiff;
          verboseStr = [ verboseStr, ',  Rel Diff: ', num2str( relDiff ) ];   %#ok<AGROW>
        end
      end

      if verbose == true, disp( verboseStr ); end
    end

    xTilde = x;
  end

  xStar = x;
end

