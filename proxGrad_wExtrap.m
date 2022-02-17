
function [xStar,objectiveValues,relDiffs,extrapolated] = proxGrad_wExtrap( x, gGrad, proxth, varargin )
  % [xStar,objectiveValues,relDiffs] = proxGrad_wExtrap( x, gGrad, proxth [, ...
  %   'g', g, 'h', h, 'N', N, 'q', q, 't', t, 'tol', tol, 'verbose', verbose ] )
  %
  % This function implements the proximal gradient optimization algorithm
  % which minimizes functions of form g(x) + h(x) where g is differentiable
  % and h has a simple proximal operator.  This algorithm includes an extrapolation
  % based on "Geometry of First Order Methods and Adaptive Acceleration" by Poon
  % and Liang
  %
  % Inputs:
  % x - the starting point
  % g - a function handle representing the g function; accepts a vector x
  %     as input and returns a scalar.
  % gGrad - a function handle representing the gradient function of g;
  %     input: the point to evaluation, output: the gradient vector
  % proxth - the proximal operator of the h function (with parameter t);
  %     two inputs: the vector and the scalar value of the parameter t
  %
  % Optional Inputs:
  % h - a handle to the h function.  This is needed to calculate the
  %     objective values.
  % N - the number of iterations that FISTA will perform (default is 100)
  % t - step size (default is 1)
  % verbose - if set then prints fista iteration
  %
  % Outputs:
  % xStar - the optimal point
  %
  % Written by Nicholas Dwork - Copyright 2017
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  defaultN = 100;
  nMaxExtrap = 50;
  defaultTol = 1d-4;

  p = inputParser;
  p.addParameter( 'g', [] );
  p.addParameter( 'h', [] );
  p.addParameter( 'N', defaultN, @(x) ispositive(x) || numel(x)==0 );
  p.addParameter( 'printEvery', 1, @ispositive );
  p.addParameter( 'q', 3, @ispositive );
  p.addParameter( 't', 1, @isnumeric );
  p.addParameter( 'tol', defaultTol, @(x) numel( x ) == 0  ||  ispositive( x ) );
  p.addParameter( 'verbose', 0, @(x) isnumeric(x) || islogical(x) );
  p.parse( varargin{:} );
  g = p.Results.g;
  h = p.Results.h;
  N = p.Results.N;  % total number of iterations
  q = p.Results.q;
  printEvery = p.Results.printEvery;  % display result printEvery iterations
  t = p.Results.t;  % t0 must be greater than 0
  tol = p.Results.tol;
  verbose = p.Results.verbose;

  if numel( N ) == 0, N = defaultN; end
  if numel( tol ) == 0,  tol = defaultTol;  end

  if t <= 0, error('fista: t0 must be greater than 0'); end

  if nargout > 1, objectiveValues = zeros(N,1); end
  if nargout > 2, relDiffs = zeros(N,1); end
  if nargout > 3, extrapolated = zeros(N,1); end

  calculateObjectiveValues = false;
  if nargout > 1
    if numel(h) == 0 || numel(g) == 0
      error( 'fista.m - Cannot calculate objective values without g and h function handles' );
    else
      calculateObjectiveValues = true;
    end
  end

  calculateRelDiffs = false;
  if nargout > 2 || ( tol > 0  &&  tol < Inf )
    calculateRelDiffs = true;
  end

  z = x;
  V = zeros( numel(z), q );
  nSinceExtrap = 0;

  k = 0;
  relDiff = Inf;
  while k < N  &&  relDiff > tol
    x = z - t * gGrad( z );
    
    lastZ = z;
    z = proxth( x, t );
    dz = z - lastZ;

    objectiveValues( k + 1 ) = g( z ) + h( z );

    if calculateRelDiffs == true
      relDiff = norm( dz(:) ) / norm( lastZ(:) );
      if nargout > 2, relDiffs( k + 1 ) = relDiff; end
    end

    if verbose>0 && mod( k, printEvery ) == 0
      formatString = ['%', num2str(ceil(log10(N))), '.', num2str(ceil(log10(N))), 'i' ];
      verboseString = [ 'proxGrad_wExtrap Iteration: ', num2str(k,formatString) ];
      if calculateObjectiveValues > 0
        verboseString = [ verboseString, ',  objective: ', num2str( objectiveValues(k+1) ) ];   %#ok<AGROW>
      end
      if calculateRelDiffs > 0
        verboseString = [ verboseString, ',  relDiff: ', num2str( relDiffs(k+1) ) ];  %#ok<AGROW>
      end
      verboseString = [ verboseString, ',  nSinceExtrap: ', num2str( nSinceExtrap ) ];  %#ok<AGROW>
      disp( verboseString );
    end

    nExtraped = 0;
    if nSinceExtrap >= q

      c = V \ dz;

      extrapRelDiff  = Inf;
      while nExtraped < nMaxExtrap  &&  extrapRelDiff  > tol  &&  k < N
        dz = V * c;
        xBar = z + dz;
        zBar = proxth( xBar, t );
        objValueBar = g( zBar ) + h( zBar );
        dObjValue = objectiveValues( k + 1 ) - objectiveValues( k );
          % Note:  k >= q at this point.
        if objValueBar - objectiveValues( k + 1 )  >  1d-6
          break;
        end

        V = circshift( V, [0 1] );
        V(:,1) = dz;

        z = zBar;
        k = k + 1;
        if nargout > 3, extrapolated( k ) = 1; end
        objectiveValues( k + 1 ) = objValueBar;
        if calculateRelDiffs == true
          extrapRelDiff = norm( dz(:) ) / norm( z(:) );
          if nargout > 2, relDiffs( k + 1 ) = extrapRelDiff; end
        end
        nExtraped = nExtraped + 1;
      end
    end

    if nExtraped == 0
      V = circshift( V, [0 1] );
      V(:,1) = dz;

      nSinceExtrap = nSinceExtrap + 1;

    else
      nSinceExtrap = 0;
    end

    k = k + 1;
  end

  if calculateObjectiveValues == true, objectiveValues = objectiveValues( 1 : k ); end
  if calculateRelDiffs == true, relDiffs = relDiffs( 1 : k ); end

  xStar = z;
end

