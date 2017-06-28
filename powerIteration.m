

function [nrm,flag] = powerIteration( M, varargin )
  % Determines an approximation for the norm of the matrix M
  %
  % [nrm,flag] = powerIteration( M, symmetric [, x0, 'maxIters', ...
  %    maxIters, 'tolerance', tolerance ])
  %
  % Inputs:
  % M is either a matrix or a function handle.
  %   If M is a function handle, then x0 must be supplied
  %   If M is a function handle, then 
  %     If M is not symmetric, it is assumed it accepts two arguments
  %       The first is the vector apply M to.
  %       The second is a string argument.  If 'transp' is supplied then the
  %         adjoint of M is appplied.
  %     If M is symmetric, then it need only accept one argument, the
  %       vector to apply M to.
  % symmetric - either true of false
  %   true if M is symmetric, false otherwise
  %
  % Optional Inputs:
  % maxIters - The maximum number of iterations permitted
  % tolerance - if the residual is less than this tolerance, then return
  %
  % Outputs:
  % nrm - approximation of norm of M
  % flag - 0 if converged; 1 if maximum iterations reached
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It 
  % is offered without any warranty expressed or implied, including the 
  % implied warranties of merchantability or fitness for a particular 
  % purpose.

  defaultMaxIters = 1000;
  defaultTolerance = 1d-4;
  p = inputParser;
  p.addRequired( 'symmetric' );
  p.addOptional( 'x0', [] );
  p.addParameter( 'maxIters', defaultMaxIters, @isnumeric );
  p.addParameter( 'tolerance', defaultTolerance, @isnumeric );
  p.parse( varargin{:} );
  symmetric = p.Results.symmetric;
  maxIters = p.Results.maxIters;
  tolerance = p.Results.tolerance;

  if isa( M, 'function_handle' )
    x0 = p.Results.x0;
    if isempty( x0 )
      error('x0 is required if M is a function handle');
    end

    if symmetric == true
      [nrm,flag] = pI_fhSymm( M, x0, maxIters, tolerance );
    else
      [nrm,flag] = pI_fh( M, x0, maxIters, tolerance );
    end

  else
    p.parse( varargin{:} );
    x0 = p.Results.x0;
    if isempty( x0 )
      x0 = rand( size(M,2), 1 );
    end
    if symmetric == true
      [nrm,flag] = pI_matSymm( M, x0, maxIters, tolerance );
    else
      [nrm,flag] = pI_mat( M, x0, maxIters, tolerance );
    end
  end
end


function [nrm,flag] = pI_fh( applyM, x, maxIters, tolerance )

  lambda = 0;
  flag = 1;
  for iter = 1:maxIters
    MtMx = applyM( applyM(x), 'transp' );
    lambdaPrev = lambda;
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


function [nrm,flag] = pI_fhSymm( applyM, x, maxIters, tolerance )

  lambda = 0;
  flag = 1;
  for iter = 1:maxIters
    Mx = applyM(x);
    lambdaPrev = lambda;
    lambda = norm(Mx,2);
    if lambda==0, break; end;

    x = Mx / lambda;

    diff = abs(lambda - lambdaPrev) / lambda;
    if diff < tolerance
      flag = 0;
      break;
    end
  end

  nrm = abs( lambda );
end


function [nrm,flag] = pI_mat( M, x, maxIters, tolerance )

  lambda = 0;
  flag = 1;
  for iter = 1:maxIters
    MtMx = M'*M*x;
    lambdaPrev = lambda;
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


function [nrm,flag] = pI_matSymm( M, x, maxIters, tolerance )

  lambda = 0;
  flag = 1;
  for iter = 1:maxIters
    Mx = M*x;
    lambdaPrev = lambda;
    lambda = norm(Mx,2);
    if lambda==0, break; end;

    x = Mx / lambda;

    diff = abs(lambda - lambdaPrev) / lambda;
    if diff < tolerance
      flag = 0;
      break;
    end
  end

  nrm = abs( lambda );
end

