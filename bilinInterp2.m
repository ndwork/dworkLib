
function out = bilinInterp2( X, Y, V, Xq, Yq, varargin )
  % out = bilinInterp2( X, Y, V, Xq, Yq [, 'op', 'notransp' ] )
  %
  % out = bilinInterp2( X, Y, V, Xq, Yq [, 'op', 'transp' ] )
  %
  % Performs a two dimensional bilinear interpolation.  Extrapolated values
  % are all set to 0.
  %
  % Inputs:
  % X - an ordered 1D sorted array of domain values
  % Y - an ordered 1D sorted array of domain values
  % V - if op is notransp, a 2D array of function values of size numel(Y) x numel( X )
  %     if op is transp, a 1D array of size numel( Xq )
  % Xq - a 1D array of query points in the domain
  % Yq - a 1D array of query points in the domain
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

  Nx = numel(X);  minX = min( X );  maxX = max( X );
  Ny = numel(Y);  minY = min( Y );  maxY = max( Y );
  subXq = Xq( Xq >= minX & Xq <= maxX & Yq >= minY & Yq <= maxY );  % NICK, YOU NEED TO ACCOUNT FOR XQ = MAXX AND YQ = MAXY
  subYq = Yq( Xq >= minX & Xq <= maxX & Yq >= minY & Yq <= maxY );
  nQ = numel( subXq );

  Xs = X(:) * ones(1,nQ);     xDiffs = bsxfun( @minus, transpose( subXq(:) ), Xs );
  Ys = Y(:) * ones(1,nQ);     yDiffs = bsxfun( @minus, transpose( subYq(:) ), Ys );

  xDiffs( xDiffs < 0 ) = Inf;
  [ lowerDiffsX, lowerIndxsX ] = min( xDiffs );  
  upperIndxsX = lowerIndxsX + 1;
  upperIndxsX( upperIndxsX > Nx ) = Nx;
  upperDiffsX = X( upperIndxsX ) - subXq(:);

  yDiffs( yDiffs < 0 ) = Inf;
  [ lowerDiffsY, lowerIndxsY ] = min( yDiffs );  
  upperIndxsY = lowerIndxsY + 1;
  upperIndxsY( upperIndxsY > Ny ) = Ny;
  upperDiffsY = Y( upperIndxsY ) - subYq(:);

  lowerDiffsX = lowerDiffsX(:);   lowerIndxsX = lowerIndxsX(:);
  lowerDiffsY = lowerDiffsY(:);   lowerIndxsY = lowerIndxsY(:);
  upperDiffsX = upperDiffsX(:);   upperIndxsX = upperIndxsX(:);
  upperDiffsY = upperDiffsY(:);   upperIndxsY = upperIndxsY(:);

  wDenomX = lowerDiffsX + upperDiffsX;
  wDenomY = lowerDiffsY + upperDiffsY;
  
  weightsLx = upperDiffsX ./ wDenomX;  weightsLx( upperDiffsX == 0 ) = 1;
  weightsLy = upperDiffsY ./ wDenomY;  weightsLy( upperDiffsY == 0 ) = 1;
  weightsUx = lowerDiffsX ./ wDenomX;  weightsUx( upperDiffsX == 0 ) = 0;
  weightsUy = lowerDiffsY ./ wDenomY;  weightsUy( upperDiffsY == 0 ) = 0;

  weightsLxLy = weightsLx .* weightsLy;
  weightsLxUy = weightsLx .* weightsUy;
  weightsUxLy = weightsUx .* weightsLy;
  weightsUxUy = weightsUx .* weightsUy;

  if strcmp( op, 'notransp' )
    % V are the values of the function at the domain indices of X
    out = zeros( numel( Xq ), 1 );
    out( Xq >= minX & Xq <= maxX & Yq >= minY & Yq <= maxY ) = ...
      V( sub2ind( [Ny,Nx], lowerIndxsY, lowerIndxsX ) ) .* weightsLxLy + ...
      V( sub2ind( [Ny,Nx], upperIndxsY, lowerIndxsX ) ) .* weightsLxUy + ...
      V( sub2ind( [Ny,Nx], lowerIndxsY, upperIndxsX ) ) .* weightsUxLy + ...
      V( sub2ind( [Ny,Nx], upperIndxsY, upperIndxsX ) ) .* weightsUxUy;

  else
    % V are the perviously queried interpolation values
    out = zeros( numel( Y ), numel( X ) );
    subV = V( Xq >= minX & Xq <= maxX & Yq >= minY & Yq <= maxY );
    subV = subV(:);
    vWeightsLxLy = subV .* weightsLxLy;
    vWeightsLxUy = subV .* weightsLxUy;
    vWeightsUxLy = subV .* weightsUxLy;
    vWeightsUxUy = subV .* weightsUxUy;

    indxsLxLy = sub2ind( [Ny,Nx], lowerIndxsY, lowerIndxsX );
    indxsLxUy = sub2ind( [Ny,Nx], upperIndxsY, lowerIndxsX );
    indxsUxLy = sub2ind( [Ny,Nx], lowerIndxsY, upperIndxsX );
    indxsUxUy = sub2ind( [Ny,Nx], upperIndxsY, upperIndxsX );

    for i = 1 : numel( lowerIndxsX )
      indxLxLy = indxsLxLy( i );  indxUxLy = indxsUxLy( i );
      indxLxUy = indxsLxUy( i );  indxUxUy = indxsUxUy( i );

      out( indxLxLy ) = out( indxLxLy ) + vWeightsLxLy( i );
      out( indxLxUy ) = out( indxLxUy ) + vWeightsLxUy( i );
      out( indxUxLy ) = out( indxUxLy ) + vWeightsUxLy( i );
      out( indxUxUy ) = out( indxUxUy ) + vWeightsUxUy( i );
    end
  end

end
