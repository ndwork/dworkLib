
function [qStar, objVal, iter] = lm( fun, q, x, y, varargin )
  % Levenberg-Marquardt algorithm for fitting f( q, x ) to ydata.  LM minimizes sum_i ( f( p, x_i ) - y_i )^2
  %
  %
  % Inputs:
  %   fun - Function handle returning f = fun( q, xdata ), where f is the model prediction
  %   p0 - Initial parameter guess (vector)
  %   xdata - Independent variable (e.g., x-coordinates or time points)
  %   ydata - Observed data (vector)
  %
  % Optional Inputs:
  %   h - step size for numerical Jacobian estimate
  %   lb - a vector of lower bounds for each parameter in q
  %   maxNumIterations - Maximum number of iterations
  %   tol - Convergence tolerance
  %   ub - a vector of upper bounds for each parameter in q
  %
  % Outputs:
  %   qStar - Optimized parameters
  %   fval - Final sum of squared residuals
  %   iter - Number of iterations performed
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
  p.addParameter( 'h', 1e-6, @ispositive );
  p.addParameter( 'lb', [] );
  p.addParameter( 'maxNumIterations', 1000, @ispositive );
  p.addParameter( 'tol', [], @(x) isempty(x) || ispositive(x) );
  p.addParameter( 'ub', [] );
  p.parse( varargin{:} );
  h = p.Results.h;
  lb = p.Results.lb;
  maxNumIterations = p.Results.maxNumIterations;
  tol = p.Results.tol;
  ub = p.Results.ub;

  if numel( lb ) > 0
    if any( q < lb )
      warning( 'Values of q are less than lb; modifying accordingly' );
      q = max( q, lb );
    end
  end

  if numel( ub ) > 0
    if any( q > ub )
      warning( 'Values of q are greater than ub; modifying accordingly' );
      q = min( q, ub );
    end
  end


  % Initialize parameters
  q = q(:); % Ensure column vector
  n = numel( q );
  lambda = 0.01; % Initial damping parameter
  
  % Evaluate initial model predictions
  [ fOut, J ] = fun( q, x );
  r = fOut - y; % Compute residuals
  objVal = sum( r.^2 ); % Objective: sum of squared residuals

  if isempty( J )
    J = zeros(length(x), n);
    for j = 1:n
      q_h = q; q_h(j) = q_h(j) + h;
      f_h = fun(q_h, x);
      J(:, j) = (f_h - fOut) / h;
    end
  end
  
  for iter = 1 : maxNumIterations

    % Compute Hessian approximation: J'*J
    H = J' * J;
    H_lm = H + lambda * diag( diag( H ) );
    
    % Compute gradient: J'*r
    g = J' * r;
    
    % Solve for step: ( J'*J + lambda*I ) * dp = -J' * r
    %dq = -H_lm \ g;
    dq = -pinv( H_lm ) * g;

    qNew = q + dq;

    if numel( lb ) > 0
      qNew = max( qNew, lb );
    end
    if numel( ub ) > 0
      qNew = min( qNew, ub );
    end

    [ fNew, J_new ] = fun(qNew, x);
    rNew = fNew - y;
    objVal_new = sum( rNew.^2 );
    
    if isempty( J_new )
      J_new = zeros( length(x), n );
      for j = 1:n
        q_h = qNew; q_h(j) = q_h(j) + h;
        f_h = fun(q_h, x);
        J_new(:, j) = ( f_h - fNew ) / h;
      end
    end

    if objVal_new < objVal
      q = qNew; r = rNew; J = J_new; objVal = objVal_new;
      lambda = lambda / 10;
    else
      lambda = lambda * 10;
    end

    if numel( tol ) > 0  &&  ( norm(dq) < tol || abs(objVal - objVal_new) < tol )
      break;
    end
  end
  
  qStar = q;
end

