
function [x,residuals] = lsqrTV( A, b, lambda, varargin )
  % x = lsqrTV_2D( A, b, lambda [, sigma, tau, 'theta', theta, ...
  %   'x0', x0, 'nIter', nIter ] );
  %
  % Solves the following regularized least squares optimization problem
  %   minimize (1/2)|| Ax - b ||_2^2 + lambda TV(x)
  % Uses Chambolle-Pock (Primal-Dual Algorithm) based on A First-Order
  % Primal-Dual Algorithm by Malitsky and Pock
  %
  % Inputs:
  % A - a matrix
  % sigma, tau - Chambolle-Pock step sizes
  % theta - the over-relaxation parameter
  % nIter - the number of iterations that CP will perform
  %
  % Outputs:
  % x - the optimal vector x
  % residuals - optionally store the residual values
  %
  % Written by Nicholas Dwork - Copyright 2017
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  defaultNIter = 1000;
  p = inputParser;
  p.addOptional( 'sigma', [] );
  p.addOptional( 'tau', [] );
  p.addParameter( 'theta', 1, @isnumeric );
  p.addParameter( 'nIter', defaultNIter, @isnumeric );
  p.addParameter( 'x0', [] );
  p.parse( varargin{:} );
  sigma = p.Results.sigma;
  tau = p.Results.tau;
  theta = p.Results.theta;
  nIter = p.Results.nIter;
  x0 = p.Results.x0;

  applyD = @(u) [ u(2:end) - u(1:end-1); 0; ];
  applyDT = @(u) [ -u(1); u(1:end-2) - u(2:end-1); u(end-1); ];

  function out = applyA( in, op )
    if isa( A, 'function_handle' )
      out = A( in, op );
    else
      if exist('op','var') && strcmp( op, 'transp' )
        out = A' * in;
      else
        out = A * in;
      end
    end
  end

  nb = numel(b);
  function out = applyK( in, op )
    if exist('op','var') && strcmp( op, 'transp' )
      in1 = in(1:nb);
      in2 = in(nb+1:end);
      out = applyA( in1, 'transp' ) + applyDT( in2 );
    else
      out1 = applyA( in, 'notransp' );
      out2 = applyD( in );
      out = [ out1(:); out2(:); ];
    end
  end

  % Initialize variables
  if isempty( x0 )
    x = applyA( b, 'transp' );
  else
    x = x0;
  end
  xBar = x;
  z1 = rand( size(b) );  % initializing z in this way makes sure it's the right shape
  z2 = applyD( x );

  % requirement for convergence guarantee: sigma * tau * norm(K)^2 <= 1
  if isempty(sigma) && isempty(tau)
    nK = powerIteration( @applyK, 0, x );
    sigma = 1/nK;
    tau = 1/nK;
  elseif isempty(sigma)
    nK = powerIteration( @applyK, 0, x );
    tau = 1 / ( sigma * nK*nK );
  elseif isempty(tau)
    nK = powerIteration( @applyK, 0, x );
    sigma = 1 / ( tau * nK*nK );
  end

  % Chambolle-Pock iteration
  if nargout > 1, residuals=zeros(nIter,1); end;
  for i=1:nIter
    tmp1 = z1 + sigma * applyA( xBar );
    tmp2 = z2 + sigma * applyD( xBar );
    z1 = ( tmp1 - sigma*b )/( 1 + sigma );
    z2 = max( min( tmp2, lambda ), -lambda );

    lastX = x;
    KTz = applyA( z1, 'transp' ) + applyDT( z2 );
    x = x - tau * KTz;  % Since G is the zero function the prox operator of G
                        % is to set x equal to the argument of the prox operator

    xBar = x + theta * ( x - lastX );

    if nargout > 1
      residuals(i) = 0.5*norm( applyK(x) - b, 2 ) + lambda * norm( x, 1 );
    end
  end

end



