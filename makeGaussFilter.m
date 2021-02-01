
function f = makeGaussFilter( N, C )
  % f = makeGaussFilter( N, C )
  % This function makes a gaussian filter of arbitrary dimensions
  % The final filter is N(1) x N(2) x ... x N(end) in size.
  % C is the covariance matrix with units in pixels.
  %   If C is a single element, the covariance matrix is assumed to be
  %     diag([ C*C C*C C*C ])
  %   If C is an array, the covariance matrix is assumed to be
  %     diag(C.*C).
  %   Otherwise, it is assumed that C is the covariance matrix

  coords = size2imgCoordinates( N );   %#ok<NASGU>
  cmd = '[ ';
  for indx = 1 : numel( N )-1
    cmd = [ cmd, 'indxs', num2str(indx), ', ' ];   %#ok<AGROW>
  end
  cmd = [ cmd, 'indxs', num2str( numel(N) ), ' ] = ndgrid( coords{:} );' ];
  eval( cmd );

  v = [];
  cmd = 'v = [ ';
  for indx = 1 : numel( N ) - 1
    cmd = [ cmd, ' indxs', num2str(indx), '(:)''; ' ];   %#ok<AGROW>
  end
  cmd = [ cmd, ' indxs', num2str(numel(N)), '(:)''; ];' ];
  eval( cmd );
  
  if ismatrix( C )
    if C <=0, error( 'Covariance matrix is not positive definite' ); end

    f = exp( -0.5 * sum( v .* ( C \ v ), 1 ) );

  else
    if min(C) <= 0, error('Covariance matrix is not positive definite'); end
    if numel(C) == 1
      C = C * ones( numel(N), 1 );
    end

    v = bsxfun( @times, v, 1./(C.*C) );
    f = exp( -0.5 * sum( v .* v, 1 ) );

  end

  f = reshape( f, N );
  f = f ./ sum( f(:) );
end
