
function [out,flag,residual] = lsqrFISTA( A, b, varargin )
  %
  % Function finds x that minimizes the L2 norm of A*x - b
  %
  % [out,flag,residual] = lsqrFISTA( A, b [, tolerance, maxIter, x0, t ]);
  %
  % Inputs:
  % A - matrix or function handle
  % b - vector
  %
  % Optional Inputs:
  % tolerance - if within this tolerance then returns
  % maxIter - the maximum number of iterations regardless of tolerance
  % x0 - initial guess of x
  % t - FISTA step size

  defaultTolerance = 1e-6;
  defaultMaxIter = 100;
  defaultX0 = [];
  defaultT = [];
  p = inputParser;
  p.addOptional( 'tolerance', defaultTolerance );
  p.addOptional( 'maxIter', defaultMaxIter );
  p.addOptional( 'x0', defaultX0 );
  p.addOptional( 't', defaultT );
  p.parse( varargin{:} );
  tolerance = p.Results.tolerance;
  maxIter = p.Results.maxIter;
  x0 = p.Results.x0;
  t = p.Results.t;

  if isa( A, 'function_handle' )
    [out,flag,residual] = lsqrFISTA_fh( A, b, tolerance, maxIter, ...
      x0, t );
  else
    [out,flag,residual] = lsqrFISTA_matrix( A, b, tolerance, maxIter, ...
      x0, t );
  end

end


function [out,flag,residual] = lsqrFISTA_fh( applyA, b, tolerance, ...
  maxIter, x0, t )

  ATb = applyA( b, 'transp' );
  bNorm = norm( b, 2 );

  if numel(x0) == 0, x0=rand(numel(ATb),1); end;

  if numel(t) == 0
    applyM = @(x) applyA( applyA( x ), 'transp' );
    normAtA = powerIteration( applyM, applyM, x0 );
    t = 1 / normAtA * 0.95;
  end

  gradG = @(Ay,ATb) applyA( Ay, 'transp' ) - ATb;

  xN1 = x0;
  xN2 = x0;
  k = 1;
  flag = 1;
  while k < maxIter
    y = xN1 + (k-2)/(k+1)*(xN1 - xN2);
    Ay = applyA(y);
    gradGy = gradG(Ay,ATb);

    x = y - t * gradGy;

    k = k + 1;
    xN2 = xN1;
    xN1 = x;

    if ( numel(tolerance) > 0 ) && ( tolerance > 0 )
      Ax = applyA( x );
      residual = norm( Ax - b, 2 ) / bNorm;
      if residual < tolerance
        flag = 0;
        break;
      end
    end
  end

  out = x;
end


function [out,flag,residual] = lsqrFISTA_matrix( A, b, tolerance, maxIter, x0, t )
  ATb = A' * b;
  bNorm = norm( b, 2 );

  if numel(x0) == 0, x0=rand(numel(ATb),1); end;

  if numel(t) == 0
    M = A' * A;
    normAtA = powerIteration( M, x0 );
    t = 1 / normAtA * 0.95;
  end

  gradG = @(Ay,ATb) A' * Ay - ATb;

  xN1 = x0;
  xN2 = x0;
  k = 1;
  flag = 1;
  while k < maxIter
    y = xN1 + (k-2)/(k+1)*(xN1 - xN2);
    Ay = A*y;
    gradGy = gradG(Ay,ATb);

    x = y - t * gradGy;

    k = k + 1;
    xN2 = xN1;
    xN1 = x;

    if ( numel(tolerance) > 0 ) && ( tolerance > 0 )
      Ax = A*x;
      residual = norm( Ax - b, 2 ) / bNorm;
      if residual < tolerance
        flag = 0;
        break;
      end
    end
  end

  out = x;
end
