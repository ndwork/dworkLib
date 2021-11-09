
function [zStar,objValues,relDiffs] = adaptiveAccelerationOptimization( z, gGrad, varargin )
  % [zStar,objValues] = adaptiveAccelerationOptimization( z, gGrad, proxth [, ...
  %   'a', a, 'alg', alg, 'b', b, 'delta', delta, 'g', g, 'h', h, 'N', N, 'q', q, ... );
  %   's', s, 't', t, 'verbose', verbose ] )
  %
  % Implements adaptive acceleration for first order methods as described in
  % "Geometry of First-Order Methods and Adaptive Acceleration" by Poon and Liang.
  % This is an implementation of algorithm 5 in the paper.
  %
  % Inputs:
  % z - initial value of optimization variable
  % gGrad - function handle to function that returns the gradient of g
  %
  % Optional Inputs:
  % alg - the first order method to acclerate.
  %   By default, if 'proxth' is supplied then proxGrad is the default algorithm,
  %   otherwise, 'gradDescent' is the default.
  %   Options are 'gradDescent', 'proxGrad', 'fista'
  %
  % Written by Nicholas Dwork, Copyright 2021
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp( 'Usage:  ' );
    disp( '  [ zStar, objValues, relDiffs ] = adaptiveAccelerationOptimization( ...' );
    disp( '    z, gGrad [, ''proxth'', proxth, ''a'', a, ''alg'', alg, ''b'', b, ...' );
    disp( '    ''delta'', delta, ''g'', g, ''h'', h, ''N'', N, ''q'', q, ''s'', s, ...' );
    disp( '    ''t'', t, ''tol'', tol, ''verbose'', verbose ] )' );
    if nargout > 0, zStar = []; end
    if nargout > 1, objValues = []; end
    if nargout > 2, relDiffs = []; end
    return;
  end

  p = inputParser;
  p.addParameter( 'a', 1, @ispositive );
  p.addParameter( 'alg', [] );
  p.addParameter( 'b', 1, @ispositive );
  p.addParameter( 'delta', 1, @ispositive );
  p.addParameter( 'g', [] );
  p.addParameter( 'h', [] );
  p.addParameter( 'N', 100, @ispositive );
  p.addParameter( 'proxth', [] );
  p.addParameter( 'q', 5, @ispositive );
  p.addParameter( 's', 10, @ispositive );
  p.addParameter( 't', 1, @ispositive );
  p.addParameter( 'tol', 1d-3, @ispositive );
  p.addParameter( 'verbose', false, @(x) numel(x)==0 || islogical(x) );
  p.parse( varargin{:} );
  a = p.Results.a;
  alg = p.Results.alg;
  b = p.Results.b;
  g = p.Results.g;
  h = p.Results.h;
  N = p.Results.N;
  proxth = p.Results.proxth;
  q = p.Results.q;
  s = p.Results.s;
  t = p.Results.t;
  tol = p.Results.tol;
  verbose = p.Results.verbose;

  if numel( alg ) == 0
    if numel( proxth ) == 0
      alg = 'gradDescent';
    else
      alg = 'proxGrad';
    end
  end

  V = zeros( numel(z), q );

  k = 0;
  nSinceLastExtrap = 0;

  if nargout > 1, objValues = zeros( N, 1 ); end

  calculateRelDiffs = false;
  if nargout > 2, calculateRelDiffs = true; end
  if numel( tol ) > 0  &&  tol > 0  &&  tol < Inf
    calculateRelDiffs = true;
  end

  if calculateRelDiffs == true, relDiffs = zeros( N, 1 ); end

  H = zeros( q );
  H( 1:q-1, 2:end ) = eye( q-1 );
  if strcmp( alg, 'fista' ), y=z; end

  relDiff = 2 * tol + 1;

  while k <= N  &&  relDiff > tol

    lastZ = z;

    %-- Update z with chosen algorithm
    if strcmp( alg, 'fista' )
      x = z - t * gGrad( x );
      lastY = y;
      y = proxth( x, t );
      z = y + ( k / (k+3) ) * ( y - lastY );

    elseif strcmp( alg, 'gradDescent' )
      z = z - t * gGrad( z );

    elseif strcmp( alg, 'proxGrad' )
      tmp = z - t * gGrad( z );
      z = proxth( tmp, t );

    else
      error( 'Unrecognized first order method' );
    end


    vNew = z - lastZ;
    zExtrapolated = false;

    if nSinceLastExtrap >= q

      c = V \ vNew;

      H(:,1) = c;

      if norm( H ) < 1
        Cs = H;
        for cPow = 2 : s
          Cs = Cs + Cs^cPow;
        end
        E = V * Cs;

        thisA = min( a, b / ( k^(1+delta) * norm( E(:) ) ) );
        z = z + thisA * E;

        zExtrapolated = true;
      end
    end

    if zExtrapolated == true
      nSinceLastExtrap = 0;
      k = k + s;

    else
      V = circshift( V, [ 0 1 ] );
      V(:,1) = vNew;

      nSinceLastExtrap = nSinceLastExtrap + 1;
      k = k + 1;
    end

    if nargout > 1
      objValues( k ) = g( z ) + h( z );
    end

    if calculateRelDiffs == true
      zDiff = z - lastZ;
      relDiff = norm( zDiff(:) ) / norm( z(:) );
      relDiffs( k ) = relDiff;
    end

    if verbose == true
      verboseStr = [ 'aao iteration: ', indx2str(k,N), ' of ', num2str(N) ];
      if nargout > 1
        verboseStr = [ verboseStr, ',  objective: ', num2str( objValues(k) ) ];   %#ok<AGROW>
      end
      if nargout > 2
        verboseStr = [ verboseStr, ',  relDiff: ', num2str( relDiffs(k) ) ];   %#ok<AGROW>
      end
      disp( verboseStr );
    end
    
  end
end


