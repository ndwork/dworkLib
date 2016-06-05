
function [nrm,flag] = powerIteration( M, varargin )
  % Determines an approximation for the norm of the matrix M
  %
  % [nrm,lambdaVals,flag] = powerIteration( M [, x0, maxIters, tolerance ])
  % [nrm,lambdaVals,flag] = powerIteration( M, MT [, x0, maxIters, tolerance ])
  %
  % Inputs:
  % M is either a matrix or a function handle.
  %   If M is a file handle, then MT (the transpose of M) must also be
  %     supplied
  %   M( x ) results in M*x;
  % 
  % Outputs:
  % nrm - approximation of norm of M
  % flag - 0 if converged; 1 if maximum iterations reached

  defaultX0 = [];
  defaultMaxIters = 500;
  defaultTolerance = 1d-4;
  p = inputParser;

  if isa( M, 'function_handle' )
    p.addRequired( 'MT' );
    p.addOptional( 'x0', defaultX0 );
    p.addOptional( 'maxIters', defaultMaxIters, @isnumeric );
    p.addOptional( 'tolerance', defaultTolerance, @isnumeric );
    p.parse( varargin{:} );
    MT = p.Results.MT;
    x0 = p.Results.x0;
    maxIters = p.Results.maxIters;
    tolerance = p.Results.tolerance;
    [nrm,flag] = powerIteration_fh( M, MT, x0, maxIters, tolerance );
  else
    p.addOptional( 'x0', defaultX0 );
    p.addOptional( 'maxIters', defaultMaxIters, @isnumeric );
    p.addOptional( 'tolerance', defaultTolerance, @isnumeric );
    p.parse( varargin{:} );
    x0 = p.Results.x0;
    maxIters = p.Results.maxIters;
    tolerance = p.Results.tolerance;
    [nrm,flag] = powerIteration_mat( M, x0, maxIters, tolerance );
  end
end


function [nrm,flag] = powerIteration_fh( ...
  applyM, applyMT, x, maxIters, tolerance )

  lambdaPrev = 0;
  flag = 1;
  for iter = 1:maxIters
    MtMx = applyMT( applyM(x) );
    lambda = norm(MtMx,2);
    if lambda==0, break; end;

    x = MtMx / lambda;

    diff = abs(lambda - lambdaPrev) / lambda;
    if diff < tolerance
      flag = 0;
      break;
    end
  end

  nrm = sqrt( lambda );
end


function [nrm,flag] = powerIteration_mat( ...
  M, x, maxIters, tolerance )

  if numel(x) == 0
    sM = size(M);
    x = rand(sM(2),1);
  end

  lambdaPrev = 0;
  flag = 1;
  for iter = 1:maxIters
    MtMx = M'*M*x;
    lambda = norm(MtMx,2);
    if lambda==0, break; end;

    x = MtMx / lambda;

    diff = abs(lambda - lambdaPrev) / lambda;
    if diff < tolerance
      flag = 0;
      break;
    end
  end

  nrm = sqrt( lambda );
end
