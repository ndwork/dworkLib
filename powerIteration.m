
function [nrm,flag] = powerIteration( M, varargin )
  % Determines an approximation for the norm of the matrix M
  %
  % [nrm,flag] = powerIteration( M [, x0, 'maxIters', maxIters, 'symmetric', true/false, ...
  %   'tolerance', tolerance, 'verbose', true/false ])
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
  %   true if M is symmetric, false (default) otherwise
  %
  % Optional Inputs:
  % x0 - initial vector for iterations
  % maxIters - The maximum number of iterations permitted (default is 1000)
  % tolerance - if the residual is less than this tolerance, then return
  % verbose - print out statements
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

  if nargin < 1
    disp( 'Usage: [nrm,flag] = powerIteration( M [, x0, ''maxIters'', maxIters, ' );
    disp( '  ''symmetric'', true/false, ''tolerance'', tolerance ])' );
    return
  end

  defaultMaxIters = 1000;
  defaultTolerance = 1d-4;
  p = inputParser;
  p.addOptional( 'x0', [] );
  p.addParameter( 'maxIters', defaultMaxIters, @isnumeric );
  p.addParameter( 'symmetric', false, @islogical );
  p.addParameter( 'tolerance', defaultTolerance, @isnumeric );
  p.addParameter( 'verbose', false, @islogical );
  p.parse( varargin{:} );
  x0 = p.Results.x0;
  maxIters = p.Results.maxIters;
  symmetric = p.Results.symmetric;
  tolerance = p.Results.tolerance;
  verbose = p.Results.verbose;

  if numel( x0 ) > 0 && max( abs( x0(:) ) ) == 0
    error( 'x0 must be nonzero' );
  end
  
  if verbose == true
    disp( 'Starting powerIteration.' );
  end

  if isa( M, 'function_handle' )
    if isempty( x0 )
      error('x0 is required if M is a function handle');
    end

    if symmetric == true
      [nrm,flag] = pI_fhSymm( M, x0, maxIters, tolerance, verbose );
    else
      [nrm,flag] = pI_fh( M, x0, maxIters, tolerance, verbose );
    end

  else
    if isempty( x0 )
      x0 = rand( size(M,2), 1 );
    end
    if symmetric == true
      [nrm,flag] = pI_matSymm( M, x0, maxIters, tolerance, verbose );
    else
      [nrm,flag] = pI_mat( M, x0, maxIters, tolerance, verbose );
    end
  end
end


function [nrm,flag] = pI_fh( applyM, x, maxIters, tolerance, verbose )
  lambda = 0;
  flag = 1;
  for iter = 1:maxIters
    MtMx = applyM( applyM(x, 'notransp'), 'transp' );
    lambdaPrev = lambda;
    lambda = norm( MtMx(:), 2 );
    if lambda==0, break; end

    if verbose == true
      disp([ 'powerIteration working on ', indx2str( iter, maxIters ), ' of ', num2str(maxIters), ...
        '   current norm estimate: ', num2str( sqrt( lambda ) ) ]);
    end

    x = MtMx / lambda;

    diff = abs(lambda - lambdaPrev) / lambda;
    if diff < tolerance
      if verbose == true, disp( 'powerIteration finished early' ); end
      flag = 0;
      break;
    end
  end

  nrm = sqrt( lambda );
end


function [nrm,flag] = pI_fhSymm( applyM, x, maxIters, tolerance, verbose )
  lambda = 0;
  flag = 1;
  for iter = 1:maxIters
    Mx = applyM(x);
    lambdaPrev = lambda;
    lambda = norm( Mx(:) );
    if lambda==0, break; end

    if verbose == true
      disp([ 'powerIteration working on ', indx2str( iter, maxIters ), ' of ', num2str(maxIters), ...
        '   current norm estimate: ', num2str( sqrt( lambda ) ) ]);
    end

    x = Mx / lambda;

    diff = abs(lambda - lambdaPrev) / lambda;
    if diff < tolerance
      if verbose == true, disp( 'powerIteration finished early' ); end
      flag = 0;
      break;
    end
  end

  nrm = abs( lambda );
end


function [nrm,flag] = pI_mat( M, x, maxIters, tolerance, verbose )
  lambda = 0;
  flag = 1;
  for iter = 1:maxIters
    MtMx = M'*M*x;
    lambdaPrev = lambda;
    lambda = norm(MtMx,2);
    if lambda==0, break; end

    if verbose == true
      disp([ 'powerIteration working on ', indx2str( iter, maxIters ), ' of ', num2str(maxIters), ...
        '   current norm estimate: ', num2str( sqrt( lambda ) ) ]);
    end

    x = MtMx / lambda;

    diff = abs(lambda - lambdaPrev) / lambda;
    if diff < tolerance
      if verbose == true, disp( 'powerIteration finished early' ); end
      flag = 0;
      break;
    end
  end

  nrm = sqrt( lambda );
end


function [nrm,flag] = pI_matSymm( M, x, maxIters, tolerance, verbose )
  lambda = 0;
  flag = 1;
  for iter = 1:maxIters
    Mx = M*x;
    lambdaPrev = lambda;
    lambda = norm(Mx,2);
    if lambda==0, break; end

    if verbose == true
      disp([ 'powerIteration working on ', indx2str( iter, maxIters ), ' of ', num2str(maxIters), ...
        '   current norm estimate: ', num2str( sqrt( lambda ) ) ]);
    end

    x = Mx / lambda;

    diff = abs(lambda - lambdaPrev) / lambda;
    if diff < tolerance
      if verbose == true, disp( 'powerIteration finished early' ); end
      flag = 0;
      break;
    end
  end

  nrm = abs( lambda );
end

