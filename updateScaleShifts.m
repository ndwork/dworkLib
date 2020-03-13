
function newParams = updateScaleShifts( g, gPrime, xs, ys, beta, alpha, ...
  domainShift, stepSize, varargin )
  % newParams = updateScaleShiftParams( ...
  %   g, gPrime, xs, ys, beta, alpha, gamma, stepSize [, 'objective', objective, ...
  %   'Lb', Lb, 'Ub', Ub, 'findBeta', findBeta, 'findAlpha', findAlpha, ...
  %   'findDomainShift', findDomainShift ] )
  %
  % Compute a projected (sub)gradient descent update on the scaling and shifting
  % parameters to minimize || ys - beta * g( alpha * xs - domainShift ) ||_2
  %
  % The parameters are beta, alpha, and gamma
  % xs - the domain values of the gamma variate function
  % ys - the noisy outputs of the shifted and scaled gamma variate function
  %
  % Inputs:
  % f - a function handle for the function itself
  % fPrime - a function handle for the derivative of the function f
  % xs - a 1D array specifying the domain variables
  % ys - a 1D array specifying the noisy ouputs
  % beta - a scalar variable
  % alpha - a scalar variable
  % domainShift - a scalar variable
  % stepSize - (sub)gradient descent step size
  %
  % Optional Inputs:
  % objective - specifies which objective function to minimize
  %   Options are: L2Sq - minimize (1/2) || beta * f( alpha t - domainShift ) - ys ||_2^2
  %                L1 - minimize | beta * f( alpha t - domainShift ) - ys |_1
  %   The default is L2Sq
  % 'Lb' - an array equal to the number of parameters specifying the lower bounds
  % 'Ub' - an array equal to the number of parameters specifying the upper bounds
  % findBeta - determine the new value of the range scaling parameter
  % findAlpha - determine the new value of the domain scaling parameter
  % findDomainShift - determine the new value of the domain shifting parameter
  %
  % Outputs:
  % newParams = an array of size 1-3 (depending on the parameters sought)
  %
  % Written by Nicholas Dwork, Copyright 2019
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addParameter( 'findBeta', 1, @(x) isnumeric(x) || islogical(x) );
  p.addParameter( 'findAlpha', 1, @(x) isnumeric(x) || islogical(x) );
  p.addParameter( 'findDomainShift', 1, @(x) isnumeric(x) || islogical(x) );
  p.addParameter( 'Lb', [], @isnumeric );
  p.addParameter( 'Ub', [], @isnumeric );
  p.addParameter( 'objective', 'L2Sq', @(x) true );
  p.parse( varargin{:} );
  findBeta = p.Results.findBeta;
  findAlpha = p.Results.findAlpha;
  findDomainShift = p.Results.findDomainShift;
  Lb = p.Results.Lb;
  Ub = p.Results.Ub;
  objective = p.Results.objective;

  if numel( stepSize ) == 1, stepSize = stepSize * ones(3,1); end

  tmp = alpha .* xs(:) - domainShift;
  gTmp = g( tmp );
  gPrimeTmp = gPrime( tmp );

  errs = beta * gTmp - ys(:);
  if strcmp( objective, 'L1' ), sErrs = sign( errs ); end

  newParams = [];
  paramIndx = 1;

  if findBeta ~= 0
    if strcmp( objective, 'L2Sq' )
      grad1  = sum( errs .* gTmp );
    elseif strcmp( objective, 'L1' )
      grad1  = sum( sErrs .* gTmp );
    end
    newBeta = beta - stepSize(paramIndx) * grad1;
    if numel( Lb ) ~= 0
      newBeta = max( newBeta, Lb(paramIndx) );
    end
    if numel( Ub ) ~= 0
      newBeta = min( newBeta, Ub(paramIndx) );
    end
    newParams = newBeta;
    paramIndx = paramIndx + 1;
  end

  if findAlpha ~= 0
    if strcmp( objective, 'L2Sq' )
      grad2 = sum( beta .* errs .* gPrimeTmp .* xs(:) );
    elseif strcmp( objective, 'L1' )
      grad2 = sum( beta .* sErrs .* gPrimeTmp .* xs(:) );
    end
    newAlpha = alpha - stepSize(paramIndx) * grad2;
    if numel( Lb ) ~= 0
      newAlpha = max( newAlpha, Lb(paramIndx) );
    end
    if numel( Ub ) ~= 0
      newAlpha = min( newAlpha, Ub(paramIndx) );
    end
    newParams = [ newParams, newAlpha ];
    paramIndx = paramIndx + 1;
  end

  if findDomainShift ~= 0
    if strcmp( objective, 'L2Sq' )
      grad3 = sum( -beta * errs .* gPrimeTmp );
    elseif strcmp( objective, 'L1' )
      grad3 = sum( -beta * sErrs .* gPrimeTmp );
    end
    newDomainShift = domainShift - stepSize(paramIndx) * grad3;
    if numel( Lb ) ~= 0
      newDomainShift = max( newDomainShift, Lb(paramIndx) );
    end
    if numel( Ub ) ~= 0
      newDomainShift = min( newDomainShift, Ub(paramIndx) );
    end
    newParams = [ newParams, newDomainShift ];
  end

end

