
function [roots,roots2] = quadRoots( p, b, c )
  % roots = quadRoots( a, b, c ) or
  % roots = quadRoots( p ) or
  % [roots1,roots2] = quadRoots( a, b, c ) or
  % [roots1,roots2] = quadRoots( p );
  %
  % Numerically stable method of Determining the roots of a quadratic polynomial
  % Written according to equation 5.6.4 of Numerical Recipes in C
  %
  % Inputs:
  %   Quadratic is either a x^2 + b x + c or
  %     p(1) + p(2) x + p(3) x^2
  % a - a 1D array (or scalar)
  % b - a 1D array (or scalar)
  % c - a 1D array (or scalar)
  % p - a 3xN array (where N is the number of quadratics to solve)
  %   Polynomials are p(1,:) + p(2,:) x + p(3,:) x^2
  %
  % Written by Nicholas Dwork, Copyright 2018
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp( 'Usage: ' );
    disp( '   roots = quadRoots( a, b, c ) or' );
    disp( '   roots = quadRoots( p ) or' );
    disp( '   [roots1,roots2] = quadRoots( a, b, c ) or' );
    disp( '   [roots1,roots2] = quadRoots( p )' );
    return
  end

  if nargin < 3
    d = p(2,:) .* p(2,:) - 4 .* p(1,:) .* p(3,:);  % discriminant
    sqrtd = sqrt( d );
    sb = sign( p(2,:) );
    sb( d<0 ) = -sign( real( conj( p(2,d<0) ) .* sqrtd( d<0 ) ) );
    sb( sb==0 ) = 1;

    q = -0.5 * ( p(2,:) + sb .* sqrtd );
    if nargout > 1
      roots = q ./ p(3,:);   roots2 = p(1,:) ./ q;
    else
      roots = [ q ./ p(3,:);  p(1,:) / q; ];
    end

  else
    % a is p
    d = b .* b - 4 .* p .* c;  % discriminant
    sqrtd = sqrt( d );
    sb = sign( b );
    sb( d<0 ) = -sign( real( conj( b( d<0 ) ) .* sqrtd( d<0 ) ) );
    sb( sb==0 ) = 1;
    q = -0.5 * ( b + sb .* sqrtd );
    
    if nargout > 1
      roots = q ./ p;   roots2 = c ./ q;
    else
      roots = [ q ./ p;  c ./ q; ];
    end
  end

end