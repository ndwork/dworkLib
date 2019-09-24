
function [x,residuals] = lassoCPwLS( K, b, gamma, varargin )
  % x = lasso_CPwLS( K, b, gamma [, 'nIter', nIter ] );
  %
  % Minimizes the Lasso problem:
  %   minimize (1/2)|| K x - b ||_2 + gamma || x ||_1
  % Uses Chambolle-Pock (Primal-Dual Algorithm) with line search based on A First-Order
  % Primal-Dual Algorithm with Linesearch by Malitsky and Pock
  %
  % Inputs:
  % nIter - the number of iterations that CP will perform
  %
  % Outputs:
  % x - the optimal x
  % res - optionally store the residual values
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  defaultNIter = 1000;
  p = inputParser;
  p.addParameter( 'nIter', defaultNIter, @isnumeric );
  p.parse( varargin{:} );
  nIter = p.Results.nIter;

  tau = sqrt(min(size(K))) / norm(K,'fro');
  mu = 0.7;
  delta = 0.99;
  beta = 1;
  theta = 1;

  x = rand( size(K,2), 1 );
  y = K*x;  % initializing z in this way makes sure it has the right shape
  if nargout > 1, residuals=zeros(nIter,1); end
  for i=1:nIter
    % step 1
    tmp = x - tau*K'*y;
    lastX = x;
    x = tmp - max( min( tmp, tau*gamma), -tau*gamma );

    % step 2
    lastTau = tau;
    tau = tau * sqrt( 1 + theta );
    lastY = y;
    dx = x - lastX;
    while true
      theta = tau / lastTau;
      xBar = x + theta * dx;

      betaTau = beta * tau;
      tmp = lastY + betaTau*K*xBar;
      y = ( tmp - betaTau * b )/( 1 + betaTau );

      dy = y - lastY;
      tmp1 = sqrt(beta)*tau*norm(K'*dy,2);
      tmp2 = delta * norm(dy,2);
      if tmp1 <= tmp2
        break;
      else
        tau = tau * mu;
      end
    end

    if nargout > 1, residuals(i) = 0.5*norm( K*x - b ) + gamma * norm( x, 1 ); end
  end

end
