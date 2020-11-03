
function [x,residuals] = lassoCP( K, b, gamma, varargin )
  % [x,residuals] = lassoCP( K, b, gamma [, sigma, tau, ...
  %   theta, 'nIter', nIter ] );
  %
  % Solves the Lasso problem:
  %   minimize (1/2)|| K x - b ||_2 + gamma || x ||_1
  % Uses Chambolle-Pock (Primal-Dual Algorithm) based on A First-Order
  %   Primal-Dual Algorithm by Malitsky and Pock
  %
  % Inputs:
  % sigma, tau - Chambolle-Pock step sizes
  % theta - the over-relaxation parameter
  % nIter - the number of iterations that CP will perform
  %
  % Outputs:
  % x - the optimal x
  % residuals - optionally store the residual values
  %
  % Written by Nicholas Dwork - Copyright 2016  
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  defaultNIter = 1000;
  p = inputParser;
  p.addOptional( 'sigma', [] );
  p.addOptional( 'tau', [] );
  p.addOptional( 'theta', 1 );
  p.addParameter( 'nIter', defaultNIter, @isnumeric );
  p.parse( varargin{:} );
  sigma = p.Results.sigma;
  tau = p.Results.tau;
  theta = p.Results.theta;
  nIter = p.Results.nIter;

  % requirement: sigma * tau * norm(K)^2 <= 1
  if isempty(sigma) && isempty(tau)
    nK = powerIteration(K, 0);
    sigma = 1/nK;
    tau = 1/nK;
  elseif isempty(sigma)
    nK = powerIteration(K, 0);
    tau = 1 / ( sigma * nK*nK );
  elseif isempty(tau)
    nK = norm(K);
    sigma = 1 / ( tau * nK*nK );
  end

  x = rand( size(K,2), 1 );
  xBar = x;
  z = K * x;  % initializing z in this way makes sure it has the right shape
  if nargout > 1, residuals=zeros(nIter,1); end
  for i=1:nIter
    tmp = z + sigma * K * xBar;
    z = ( tmp - sigma * b )/( 1 + sigma );

    tmp = x - tau * K' * z;
    lastX = x;
    x = tmp - max( min( tmp, tau*gamma), -tau*gamma );

    xBar = x + theta * ( x - lastX );

    if nargout > 1
      residuals(i) = 0.5*norm( K*x - b ) + gamma * norm( x, 1 );
    end
  end

end
