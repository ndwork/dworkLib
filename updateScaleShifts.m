
function newParams = updateScaleShifts( g, gPrime, ts, ys, beta, alpha, ...
  domainShift, stepSize, varargin )
  % newParams = updateScaleShiftParams( ...
  %   g, gPrime, ts, ys, beta, alpha, gamma, stepSize [, ...
  %   'findBeta', findBeta, 'findAlpha', findAlpha, 'findGamma', findGamma ] )
  %
  % Compute a gradient descent update on the scaling and shifting parameters to minimize
  %   (1/2) || beta * f( alpha t - domainShift ) - ys ||_2^2
  % The parameters are beta, alpha, and gamma
  % ts - the domain values of the gamma variate function
  % ys - the noisy outputs of the shifted and scaled gamma variate function
  %
  % Inputs:
  % f - a function handle for the function itself
  % fPrime - a function handle for the derivative of the function f
  % ts - a 1D array specifying the domain variables
  % ys - a 1D array specifying the noisy ouputs
  % beta - a scalar variable
  % alpha - a scalar variable
  % domainShift - a scalar variable
  % stepSize - gradient descent step size
  %
  % Optional Inputs:
  % findBeta - determine the new value of the range scaling parameter
  % findAlpha - determine the new value of the domain scaling parameter
  % findGamma - determine the new value of the domain shifting parameter
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
  p.parse( varargin{:} );
  findBeta = p.Results.findBeta;
  findAlpha = p.Results.findAlpha;
  findDomainShift = p.Results.findDomainShift;

  if numel( stepSize ) == 1
    stepSize = stepSize * ones(3,1);
  end

  tmp = alpha .* ts(:) - domainShift;
  gTmp = g( tmp );
  gPrimeTmp = gPrime( tmp );

  errs = beta * gTmp - ys(:);

  newParams = [];

  if findBeta ~= 0
    grad1  = sum( errs .* gTmp );
    newBeta = beta - stepSize(1) * grad1;
    newParams = newBeta;
  end

  if findAlpha ~= 0
    grad2 = sum( beta .* errs .* gPrimeTmp .* ts(:) );
    newAlpha = alpha - stepSize(2) * grad2;
    newParams = [ newParams, newAlpha ];
  end

  if findDomainShift ~= 0
    grad3 = sum( -errs .* gPrimeTmp );
    newDomainShift = domainShift - stepSize(3) * grad3;
    newParams = [ newParams, newDomainShift ];
  end

end

