
function [p,dp] = evaluatePoly( c, x )
  % [p,dp] = evauatePoly( c, x )
  % Evaluates the polynomial defined by the coefficients in c at the
  % values of x.  This function was written according to section 5.3
  % of Numerical Recipes in C
  % 
  % Inputs:
  % c - the polynomial coefficients
  %   p(x) = c(1) + c(2)*x + c(3)*x^3 + ... + c(N)*x^(N)
  % x - a 1D array of domain variable values
  %
  % Ouptuts:
  % p - the polynomial's values
  %
  % Optional Ouptuts:
  % dp - the derivative of the polynomial
  %
  % Written by Nicholas Dwork, Copyright 2018
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = c(end);

  if nargout > 1

    % evaluate the polynomial and its derivatve
    dp = 0;
    for i = numel(c)-1:-1:1
      dp = p + dp.*x;
      p = p.*x + c(i);
    end

  else

    for i = numel(c)-1:-1:1
      p = p.*x + c(i);
    end

  end
end

