
function [x,residuals] = lsqrTV_2D( applyA, b, lambda, varargin )
  % x = lsqrTV_2D( A, b, lambda [, sigma, tau, 'theta', theta, ...
  %   'x0', x0, 'nIter', nIter ] );
  %
  % Solves the following regularized least squares optimization problem
  %   minimize (1/2)|| Ax - b ||_2^2 + lambda TV(x)
  % Uses Chambolle-Pock (Primal-Dual Algorithm) based on A First-Order
  % Primal-Dual Algorithm by Malitsky and Pock
  %
  % Inputs:
  % applyA - a function handle for a function that accepts two arguments
  %   the first argument is the 2D array to be operated on
  %   the second argument is a string.  If 'transp' then the function
  %     returns the transpose of the operation.
  % b - a 2D array
  % lambda - the regularization parameter
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

  % Initialize variables
  if isempty( x0 )
    x = applyA( b, 'transp' );
  else
    x = x0;
  end
  xBar = x;
  
  [nRows,nCols] = size(x);
  applyD1 = @(u) cat( 2, u(:,2:end) - u(:,1:end-1,:), zeros(nRows,1) );
  applyD2 = @(u) cat( 1, u(2:end,:) - u(1:end-1,:), zeros(1,nCols) );
  applyD1T = @(u) cat( 2, -u(:,1), u(:,1:end-2) - u(:,2:end-1), u(:,end-1) );
  applyD2T = @(u) cat( 1, -u(1,:), u(1:end-2,:) - u(2:end-1,:), u(end-1,:) );

  z1 = rand( size(b) );  % initializing z in this way makes sure it's the right shape
  z2 = applyD1( x );
  z3 = applyD2( x );

  function out = applyK( in, op )
    if exist('op','var') && strcmp( op, 'transp' )
      in1 = in(1:size(b,1),:);
      in2 = in(size(b,1)+1:size(b,1)+size(x,1),:);
      in3 = in(size(b,1)+size(x,1)+1:end,:);
      tmp1 = applyA( in1, 'transp' );
      tmp2 = applyD1T( in2 );
      tmp3 = applyD2T( in3 );
      out = tmp1 + tmp2 + tmp3;
    else
      out1 = applyA( in, 'notransp' );
      out2 = applyD1( in );
      out3 = applyD2( in );
      out = [ out1; out2; out3; ];
    end
  end

  % requirement for convergence guarantee: sigma * tau * norm(K)^2 <= 1
nK=0;  load( 'nK.mat' );
  if isempty(sigma) && isempty(tau)
    %nK = powerIteration( @applyK, 0, x(:) );
    sigma = 1/nK;
    tau = 1/nK;
  elseif isempty(sigma)
    %nK = powerIteration( @applyK, 0, x );
    tau = 1 / ( sigma * nK*nK );
  elseif isempty(tau)
    %nK = powerIteration( @applyK, 0, x );
    sigma = 1 / ( tau * nK*nK );
  end

  % Chambolle-Pock iteration
  if nargout > 1, residuals=zeros(nIter,1); end;
  for i=1:nIter
    tmp1 = z1 + sigma * applyA( xBar );
    tmp2 = z2 + sigma * applyD1( xBar );
    tmp3 = z3 + sigma * applyD2( xBar );
    z1 = ( tmp1 - sigma*b )/( 1 + sigma );
    z2 = max( min( tmp2, lambda ), -lambda );
    z3 = max( min( tmp3, lambda ), -lambda );

    lastX = x;
    KTz = applyA( z1, 'transp' ) + applyD1T( z2 ) + applyD2T( z3 );
    x = x - tau * KTz;  % Since G is the zero function the prox operator of G
                        % is to set x equal to the argument of the prox operator

    xBar = x + theta * ( x - lastX );

    if nargout > 1
      Ax = applyA(x);
      D1x = applyD1(x);
      D2x = applyD2(x);
      residuals(i) = 0.5*norm( Ax(:) - b(:), 2 )^2 + ...
        lambda * norm( D1x(:), 1 ) + lambda * norm( D2x(:), 1 );
    end
  end

end



