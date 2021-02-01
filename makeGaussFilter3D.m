
function f = makeGaussFilter3D( N, C )
  % f = gaussFilter3D( N, C )
  % This function makes a 3D gaussian filter
  % The final filter is N(1)xN(2)xN(3) in size.  If N is a single
  %   element, then the final filters is NxNxN in size.
  % C is the covariance matrix.
  %   If C is a single element, the covariance matrix is assumed to be
  %     diag([ C*C C*C C*C ])
  %   If C is a 3 element array, the covariance matrix is assumed to be
  %     diag(C.*C).
  %   Otherwise, it is assumed that C is the covariance matrix

  error('This code looks buggy.  Either check it, or call makeGaussFilter instead.');

  if numel(N) == 1
    N = [N N N];
  end

  [y,x,z] = ind2sub( N, 1:prod(N) );
  y = y - (N(1)+1)/2;   y = reshape(y,N);
  x = x - (N(2)+1)/2;   x = reshape(x,N);
  z = z - (N(3)+1)/2;   z = reshape(z,N);

  if numel(C) == 1

    if C <= 0
      error('Covariance matrix is not positive definite');
    end
    f = exp( -0.5/(C*C) * ( x.*x + y.*y + z.*z ) );

  elseif numel(C) == 3

    if min(C) <= 0
      error('Covariance matrix is not positive definite');
    end
    f = exp( -0.5 * ( x.*x/C(1) + y.*y/C(2) + z.*z/C(3) ) );

  else

    v = [ x; y; z; ];
    Cinv = inv(C);
    f = zeros( size( v, 1 ), 1 );
    for i= 1 : numel( f )
      f(i) = exp( -0.5 * v(i,:)' * Cinv * v(i,:) );   %#ok<MINV>
    end
    f = reshape( f, N );

  end

  f = f ./ sum( f(:) );
end
