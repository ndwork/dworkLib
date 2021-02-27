
function out = linInterp( X, V, Xq, varargin )
  % out = linInterp( X, V, Xq [, 'op', 'notransp' ] )
  %
  % out = linInterp( X, V, Xq [, 'op', 'transp' ] )
  %
  % Performs a 1 dimensional linear interpolation.  Extrapolated values
  % are all set to 0.
  %
  % Inputs:
  % X - an ordered 1D array of domain values
  % V - if op is notransp, a 1D array of function values of size numel( X )
  %     if op is transp, a 1D array of size numel( Xq )
  % Xq - a 1D array of query points in the domain
  % 
  % Optional Inputs:
  % op - by default, performs linear interpolation.  If set to 'transp', performs adjoint.
  %
  % Written by Nicholas Dwork - Copyright 2021
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addParameter( 'op', 'notransp', @(x) true );
  p.parse( varargin{:} );
  op = p.Results.op;

  Nx = numel( X );  minX = min( X );  maxX = max( X );
  subXq = Xq( Xq >= minX & Xq <= maxX );
  nQ = numel( subXq );

  Xs = X(:) * ones(1,nQ);
  xDiffs = bsxfun( @minus, transpose( subXq(:) ), Xs );

  xDiffs( xDiffs < 0 ) = Inf;
  [ lowerXDiffs, lowerIndxs ] = min( xDiffs );  
  upperIndxs = lowerIndxs + 1;
  upperIndxs( upperIndxs > Nx ) = Nx;
  upperXDiffs = X( upperIndxs ) - subXq(:);

  lowerXDiffs = lowerXDiffs(:);   lowerIndxs = lowerIndxs(:);
  upperXDiffs = upperXDiffs(:);   upperIndxs = upperIndxs(:);

  wDenom = lowerXDiffs + upperXDiffs;
  lowerWeights = upperXDiffs ./ wDenom;  lowerWeights( upperXDiffs == 0 ) = 1;
  upperWeights = lowerXDiffs ./ wDenom;  upperWeights( upperXDiffs == 0 ) = 0;

  if strcmp( op, 'notransp' )
    % V are the values of the function at the domain indices of X
    out = zeros( numel( Xq ), 1 );
    out( Xq >= minX & Xq <= maxX ) = ...
      V( lowerIndxs ) .* lowerWeights + V( upperIndxs ) .* upperWeights;

  else
    % V are the perviously queried interpolation values
    out = zeros( numel( X ), 1 );
    vLowerWeights = V( Xq >= minX & Xq <= maxX ) .* lowerWeights;
    vUpperWeights = V( Xq >= minX & Xq <= maxX ) .* upperWeights;
    for i = 1 : numel( lowerIndxs )
      indxLower = lowerIndxs( i );
      indxUpper = upperIndxs( i );
      out( indxLower ) = out( indxLower ) + vLowerWeights( i );
      out( indxUpper ) = out( indxUpper ) + vUpperWeights( i );
    end
  end

end
