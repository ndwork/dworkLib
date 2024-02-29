function [xstar, iters, alphas] = dr_wLS(x0, proxf, proxg, maxIter)
  %dr_wLS Douglas-Rachford with line-search
  %
  %   Douglas-Rachford solves problems of the form
  %       min  f(x) + g(x)
  %        x
  % and is especially useful if f and g have simple-to-evaluate proximal
  % operators. This implementation of DR turns it into averaged operator
  % iteration. Given the proximal operators for f and g, the reflected
  % proximal operators are
  %
  %     Rf(x, gamma) = 2 * proxf(x, gamma) - x
  %     Rg(x, gamma) = 2 * proxg(x, gamma) - x
  %
  % and operator S is defined as
  %     Sx = Rg( Rf( x, gamma ), gamma )
  %
  % the DR iteration is then written as
  %     z_k+1 = z_k + alpha( Sz_k - z_k )
  % where alpha, the relaxation parameter, is chosen via linesearch.
  %
  % The algorithm is written as
  %     r_k    = Sx_k - x_k
  %     xbar_k = x_k + alpha_bar * r_k
  %     rbar_k = Sxbar_k - xbar_k
  %     x_k+1  = x_k + alpha_k * r_k
  %
  % The linesearch chooses alpha_k such that
  %     || r_k+1 ||_2 = || Sx_k+1 - x_k+1 ||_2 <= (1 - eps) || rbar_k ||_2
  %
  % The default behavior of this function is to do a backtracking linesearch
  % starting at alpha_max = 50 and searching using alpha_k+1 = (1 / 1.4) alpha_k
  % where these values are taken from the Boyd paper.
  % We set alpha_bar = 1/2, and allow alpha_k to be in [alpha_bar, alpha_max]
  % i.e. if the linesearch gets below alpha_bar, alpha_k = alpha_bar
  % this ensures the linesearch takes a maximum of 14 steps
  % default value of gamma is 3
  % 
  %    INPUTS
  %     x0 - starting point
  %     proxf - proximal operator of f, should accept a scaling parameter
  %             i.e. y = proxf(x, lambda)
  %     proxg - proximal operator of g, should accept a scaling parameter
  %             i.e. y = proxf(x, lambda)
  %     maxIter - maximum iterations to run
  %
  %    OUTPUTS
  %     xstar - x*, optimal point
  %     iters - || Sxk - xk ||_2 at each iteration, where S is defined above
  %     alphas - step size alphak chosen via linesearch
  %
  % Written by Alex McManus - Copyright 2024
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.


Rg = @(x, gamma) 2 * proxg(x, gamma) - x;
Rf = @(x, gamma) 2 * proxf(x, gamma) - x;

%%% parameters
alpha_bar = 0.5; % alpha_bar
gamma = 3; % gamma for prox operators
eps = 0.03; % eps for (1 - eps) || rbar_k || in linesearch
tol = 1e-3; % tolerance for exit criterion
alpha_change = 1/1.4; % factor for change in alpha during linesearch

S = @(x) Rg( Rf( x, gamma ), gamma );

iters = zeros(maxIter, 1);
alphas = zeros(maxIter, 1);

xk = x0;

for i = 1:maxIter
%   if mod(i, 10) == 0
%     fprintf('iter %d\n', i)
%   end
  rk = S(xk) - xk;
  xk_bar = xk + alpha_bar*rk;
  rk_bar = S(xk_bar) - xk_bar;

  iters(i) = norm(rk);
  if norm(rk) < tol
    break
  end

  alpha_k = 50;
  subiter = 0;
  while true
    subiter = subiter + 1;
    %fprintf('subiter %d\n', subiter);
    xkp1 = xk + alpha_k*rk;
    rkp1 = S(xkp1) - xkp1;
    if norm(rkp1) < (1-eps)*norm(rk_bar)
      xk = xkp1;
      alphas(i) = alpha_k;
      break
    end
    alpha_k = alpha_k*alpha_change;
    if alpha_k < alpha_bar
      alpha_k = alpha_bar;
      xkp1 = xk + alpha_k*rk;
      xk = xkp1;
      alphas(i) = alpha_k;
      break
    end
  end

  xstar = xk;


end
