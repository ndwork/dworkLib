
function [xStar,oValues] = lsqrTV( applyA, b, x0, lambda, varargin )
  % xStar = lsqrTV( A, b, x0, lambda [, sigma, tau, 'theta', theta, ...
  %   'nIter', nIter ] );
  %
  % Solves the following regularized least squares optimization problem
  %   minimize (1/2)|| A x - b ||_2^2 + lambda TV(x)
  % Uses Chambolle-Pock (Primal-Dual Algorithm) based on A First-Order
  % Primal-Dual Algorithm by Malitsky and Pock
  %
  % Inputs:
  % applyA - a function handle for a function that accepts two arguments
  %   the first argument is the array to be operated on
  %   the second argument is a string.  If 'transp' then the function
  %     returns the transpose of the operation.
  % b - a 2D array
  % lambda - the regularization parameter
  % sigma, tau - Chambolle-Pock step sizes
  % theta - the over-relaxation parameter
  % nIter - the number of iterations that CP will perform
  %
  % Outputs:
  % xStar - the optimal vector x
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
  p.addParameter( 'doCheckAdjoint', false, @islogical );
  p.addParameter( 'nIter', defaultNIter, @isnumeric );
  p.addParameter( 'verbose', false, @islogical );
  p.parse( varargin{:} );
  sigma = p.Results.sigma;
  tau = p.Results.tau;
  doCheckAdjoint = p.Results.doCheckAdjoint;
  nIter = p.Results.nIter;
  verbose = p.Results.verbose;

  nb = numel( b );

  % K is the concatenation of A and the discrete gradient operator
  function out = applyK( in, op )
    if nargin < 2 || strcmp( op, 'notransp' )
      in = reshape( in, size( x0 ) );
      Ain = applyA( in );
      Din = computeGradient( in );
      out = [ Ain(:); Din(:); ];

    elseif strcmp( op, 'transp' )
      in1 = reshape( in( 1 : nb ), size(x0) );
      ATin1 = applyA( in1, 'transp' );

      in2 = reshape( in( nb + 1 : end ), [ size(x0) 2 ] );
      DTin2 = computeGradient( in2, 'transp' );

      out = ATin1 + DTin2;

    else
      error( 'Unrecognized operator for K' );
    end
  end

  if doCheckAdjoint == true
    [checkAdjK,checkErr] = checkAdjoint( x0, @applyK );
    if checkAdjK == false
      error([ 'Adjoint of K failed with error ', num2str( checkErr ) ]);
    end
  end

  f = @(x) 0;
  proxf = @(x,t) x;
  
  g1 = @(x) norm( x - b(:) )^2;
  g2 = @(x) normL2L1( x );
  g = @(x) g1( x(1:nb) ) + lambda * g2( x(nb+1:end) );

  proxg1Conj = @(x,t) proxConjL2Sq( x, t, 1, b(:) );
  proxg2Conj = @(x,t) proxConjL2L1( x, t, lambda );

  function out = proxgConj( in, t )
    out1 = proxg1Conj( in(1:nb), t );
    out2 = proxg2Conj( reshape( in(nb+1:end), [size(x0) 2] ), t );
    out = [ out1(:); out2(:); ];
  end

  useLineSearch = true;
  if useLineSearch == true
    [ xStar, oValues ] = pdhgWLS( b, proxf, @proxgConj, 'A', @applyK, ...
      'N', nIter, 'f', f, 'g', g, 'tau', tau, 'verbose', verbose );

  else  
    % requirement for convergence guarantee: sigma * tau * norm(K)^2 <= 1
    if isempty( sigma ) || isempty( tau )
      nK = powerIteration( @applyK, 0, x0(:) );

      if isempty(sigma) && isempty(tau)
        tau = 1 / nK;
      elseif isempty(sigma)
        sigma = 1 / ( tau * nK * nK );
      elseif isempty(tau)
        tau = 1 / ( sigma * nK*nK );
      end
    end

    [ xStar, oValues ] = pdhg( b, proxf, @proxgConj, tau, 'sigma', sigma, ...
      'A', @applyK, 'f', f, 'g', g, 'N', nIter, 'normA', nK, 'verbose', verbose );
    xStar = reshape( xStar, size( b ) );
    
  end

end


