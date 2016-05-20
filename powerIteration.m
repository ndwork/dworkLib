
function [nrm,lambdaVals,flag] = powerIteration( M, varargin )
  % Determines an approximation for the norm of the matrix M
  %
  % [nrm,lambdaVals,flag] = powerIteration( M [, x0, maxIters, tolerance ])
  %
  % Inputs:
  % M is either a matrix or a function handle.  M( in )
  %   in is the input vector.
  % 
  % Outputs:
  % nrm - approximation of norm of M
  % lambdaVals - 
  % flag - 0 if converged; 1 if maximum iterations reached

  defaultX0 = [];
  defaultMaxIters = 500;
  defaultTolerance = 1d-4;
  p = inputParser;
  p.addOptional( 'x0', defaultX0 );
  p.addOptional( 'maxIters', defaultMaxIters, @isnumeric );
  p.addOptional( 'tolerance', defaultTolerance, @isnumeric );
  p.parse( varargin{:} );
  x0 = p.Results.x0;
  maxIters = p.Results.maxIters;
  tolerance = p.Results.tolerance;

  if isa( M, 'function_handle' )
    [nrm,lambdaVals,flag] = powerIteration_fh( M, x0, maxIters, tolerance );
  else
    [nrm,lambdaVals,flag] = powerIteration_mat( M, x0, maxIters, tolerance );
  end
end

function [nrm,lambdaVals,flag] = powerIteration_fh( ...
  applyM, x, maxIters, tolerance )

  if numel(x)==0,
    error('Must supply an initial x0 vector if M is a file handle.');
  end

  lambdaPrev = 0;
  lambdaVals = zeros( maxIters, 1 );
  flag=1;
  for iter = 1:maxIters
     Mx = applyM( x );
     lambda = norm(Mx(:),2);
     if lambda == 0, break; end;
     x = Mx/lambda;
     lambdaVals(iter) = lambda;

     diff = abs(lambda - lambdaPrev);
     if diff < tolerance
       flag = 0;
       break;
     end

     lambdaPrev = lambda;
  end

  lambdaVals = lambdaVals(1:iter);
  nrm = sqrt(lambda);
end


function [nrm,lambdaVals,flag] = powerIteration_mat( M, x, ...
  maxIters, tolerance )

  if numel(x) == 0
    sM = size(M);
    x = rand(sM(2),1);
  end

  lambdaPrev = 0;
  lambdaVals = zeros( maxIters, 1 );
  flag = 1;
  for iter = 1:maxIters
     Mx = M*x;
     lambda = norm(Mx(:),2);
     if lambda == 0, break; end;

     x = Mx/lambda;
     lambdaVals(iter) = lambda;

     diff = abs(lambda - lambdaPrev)/lambda;
     if diff < tolerance
       flag = 0;
       break;
     end;

     lambdaPrev = lambda;
  end

  lambdaVals = lambdaVals(1:iter);
  nrm = sqrt(lambda);
end
