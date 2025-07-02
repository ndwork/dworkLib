
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

  function [ f, J ] = powerLaw( p, x )
    a = p(1); b = p(2);
    f = a * x .^ b;
    J = [x .^ b, a * x .^ b .* log(x)];   % Jacobian: [ df/da, df/db ]
  end

  % Finds a and by by fitting a line to log(y) v log(x)
  logX = log( x( x~=0 ) );
  logY = log( y( x~=0 ) );
  X = [ ones( numel(logX), 1 ), logX ];
  coeffs = X \ logY;
  a = exp( coeffs(1) );
  b = coeffs(2);

  if p.Results.linear == false
    % Uses the result of the linear fit as an initial guess.
    params = lm( @powerLaw, [a, b], x, y, 'lb', lb, 'ub', ub );
    a = params(1);
    b = params(2);
  end

  % figure;  plotnice( x, y );  hold all;  plotnice( x, a*x.^b );
  out = [ a b ];
end

