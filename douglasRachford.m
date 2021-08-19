
function out = douglasRachford( x0, proxf, proxg, t, varargin )
  % out = douglasRachford( x0, proxf, proxg, t [, 
  %   'N', N, 'rho', rho, 'verbose', verbose ] )
  %
  % minimizes f( x ) + g( x )
  %
  % Inputs:
  % x0 - an array specifying the initial input
  % proxf - a function handle to the proximal operator of f
  % proxg - a function handle to the proximal operator of g
  % t - the step size
  %
  % Optional Inputs:
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
    disp( '  out = douglasRachford( x0, proxf, proxg, t [, ' );
    disp( '    ''N'', N, ''rho'', rho, ''verbose'', verbose ] ) ' );
    if nargout > 0, out = []; end
    return
  end

  p = inputParser;
  p.addRequired( 'x0', @isnumeric );
  p.addParameter( 'N', 100, @ispositive );
  p.addParameter( 'rho', 1, @(x) numel(x) == 1 && x > 0 && x < 2 );
  p.addParameter( 'verbose', false, @islogical );
  p.parse( x0, varargin{:} );
  N = p.Results.N;
  rho = p.Results.rho;
  verbose = p.Results.verbose;

  z = x0;

  for optIter = 1 : N
    if verbose == true
      disp([ 'douglasRachford: Working on ', indx2str(optIter,N), ' of ', num2str(N) ]);
    end

    x = proxf( z, t );

    y = proxg( 2 * x - z, t );

    z = z + rho( y - x );

  end

  out = y;
end
