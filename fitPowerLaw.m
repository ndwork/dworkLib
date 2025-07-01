
function out = fitPowerLaw( x, y, varargin )
  % out = fitPowerLaw( x, y [, 'lb', lb', 'ub', ub, 'linear', true/false ] )
  %
  % Finds a and b such that y ~ a * x^b
  %
  % Optional Inputs:
  % lb - a two element array specifying the lower bounds of a and b
  % linear - if set to true, then fits a line to log(a) + b log(x).  Otherwise, solves a least-squares problem.
  % ub - a two element array specifying the upper bounds of a and b
  %
  % Outputs:
  % out = [ a b ]
  %
  % Written by Nicholas Dwork, Copyright 2025
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addParameter( 'lb', [] );
  p.addParameter( 'linear', false );
  p.addParameter( 'ub', [] );
  p.addParameter( 'verbose', false );
  p.parse( varargin{:} );
  lb = p.Results.lb;
  ub = p.Results.ub;

  % Finds a and by by fitting a line to log(y) v log(x)
  logX = log( x( x~=0 ) );
  logY = log( y( x~=0 ) );
  X = [ ones( numel(logX), 1 ), logX ];
  coeffs = X \ logY;
  a = exp( coeffs(1) );
  b = coeffs(2);

  if p.Results.linear == false

    powerLaw = @(params, x) params(1) * x .^ params(2); % params = [a, b]
    initialParams = [a, b];
    if p.Results.verbose == true
      displayOption = 'iter';
    else
      displayOption = 'none';
    end
    options = optimoptions('lsqcurvefit', ...
      'Algorithm', 'levenberg-marquardt', ...
      'Display', displayOption, ...
      'MaxIterations', 1000, ...
      'FunctionTolerance', 1e-6 );

    [params, resnorm, residual, exitflag, output] = lsqcurvefit( ...
      powerLaw, initialParams, x, y, lb, ub, options );   %#ok<ASGLU>

    a = params(1);
    b = params(2);
  end

  % figure;  plotnice( x, y );  hold all;  plotnice( x, a*x.^b );
  out = [ a b ];
end

