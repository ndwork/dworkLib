
function [ out, objValues ] = primalDualDR( x0, proxf, proxgConj, t, varargin )
  % [ out, objValues ] = primalDualDR( x0, proxf, proxgConj, t [, 
  %   'N', N, 'rho', rho, 'verbose', verbose ] )
  %
  % minimizes f( x ) + g( x )
  %
  % Inputs:
  % x0 - an array specifying the initial input
  % proxf - a function handle to the proximal operator of f
  % proxgConj - a function handle to the proximal operator of the conjugate function of g
  % t - the step size
  %
  % Optional Inputs:
  % f - the function handle for f
  % g - the function handle for g
  % N - the number of iterations to run primal-dual DR
  % rho - the relaxation parameter; 0 < rho < 2
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp( 'Usage:  ' );
    disp( '  out = primalDualDR( x0, proxf, proxgConj, t [, ' );
    disp( '    ''N'', N, ''rho'', rho, ''verbose'', verbose ] ) ' );
    if nargout > 0, out = []; end
    return
  end

  p = inputParser;
  p.addRequired( 'x0', @isnumeric );
  p.addParameter( 'f', [] );
  p.addParameter( 'g', [] );
  p.addParameter( 'N', 100, @ispositive );
  p.addParameter( 'rho', 1, @(x) numel(x) == 1 && x > 0 && x < 2 );
  p.addParameter( 'verbose', false, @islogical );
  p.parse( x0, varargin{:} );
  f = p.Results.f;
  g = p.Results.g;
  N = p.Results.N;
  rho = p.Results.rho;
  verbose = p.Results.verbose;

  if nargout > 1
    objValues = zeros( N, 1 );
  end

  z = x0;

  for optIter = 1 : N
    x = proxf( z, t );

    y = t * proxgConj( (1/t) * ( 2 * x - z ), 1/t );

    z = z + rho * ( x - z - y );

    
    if nargout > 1 || ( numel(f) > 0 && numel(g) > 0 && verbose == true )
      objValue = f( x ) + g( x );
    end

    if nargout > 1
      objValues( optIter ) = objValue;
    end

    if verbose == true
      outStr = [ 'primal-dual DR: Completed ', indx2str(optIter,N), ' of ', num2str(N) ];
      if numel( f ) > 0 && numel( g ) > 0
        outStr = [ outStr, '  objective: ', num2str( objValue ) ];   %#ok<AGROW>
      end
      disp( outStr );
    end
  end

  out = z;
end
