
function roots = quarticRoots( poly, a, b, c, d )
  % roots = quarticRoots( a, b, c, d, e ) or
  % roots = quarticRoots( b, c, d, e ) or
  % roots = quarticRoots( p )
  %
  % Solves problems of the form:  a x^4 + b x^3 + c x^2 + d x + e = 0
  % Based on "An analytic solution to Wahba's problem" by Yang and Zhou 
  %
  % Inputs:
  %   Quartic is either a x^4 + b x^3 + c x^2 + d x + e or
  %     x^4 + b x^3 + c x^2 + d x + e or 
  %     p(5,:) x^4 + p(4,:) x^3 + p(3,:) x^2 + p(2,:) x + p(1,:) or
  %     x^4 + p(4,:) x^3 + p(3,:) x^2 + p(2,:) x + p(1,:)
  %
  % Written by Nicholas Dwork, Copyright 2021
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin > 4
    a = a ./ poly;
    b = b ./ poly;
    c = c ./ poly;
    d = d ./ poly;

  elseif nargin > 3
    d = c;
    c = b;
    b = a;
    a = poly;
    clear poly;

  else
    if size( poly,1 ) == 5
      d = poly(1,:) ./ poly(5,:);
      c = poly(2,:) ./ poly(5,:);
      b = poly(3,:) ./ poly(5,:);
      a = poly(4,:) ./ poly(5,:);
    else
      d = poly(1,:);
      c = poly(2,:);
      b = poly(3,:);
      a = poly(4,:);
    end
  end

  p = a .* c - b .* b / 3 - 4 * d;
  q = a .* b .* c / 3 - a .* a .* d - 2 * b .*b .* b / 27 - c .* c + 8 * b .* d / 3;
  

  sqrtD = sqrt( (0.5*q).^2 + (p/3).^3 );  % square root of the discriminant
  pTmp = ( -0.5 * q + sqrtD ).^(1/3);
  nTmp = ( -0.5 * q - sqrtD ).^(1/3);

  y1 = pTmp + nTmp;

  % Swap so that y1 is real
  w1 = -0.5 + 0.5i * sqrt(3);
  w2 = -0.5 - 0.5i * sqrt(3);
  if sum( imag(y1(:)) > 0 ) > 0
    y2 = w1 .* pTmp + w2 .* nTmp;
    y1( imag(y1) > 0  &  imag(y2) == 0 ) = y2( imag(y1) > 0  &  imag(y2) == 0 );
  end
  if sum( imag(y1(:)) > 0 ) > 0
    y3 = w2 .* pTmp + w1 .* nTmp;
    y1( imag(y1) > 0  &  imag(y3) == 0 ) = y3( imag(y1) > 0  &  imag(y3) == 0 );
  end

  [g1,g2] = quadRoots( 1, -a, 2*b/3 - y1 );
  [h1,h2] = quadRoots( 1, -y1 - b/3, d );

  hTmp = abs( g1 .* h2 + g2 .* h1 - c ) > 0;
  if sum( hTmp(:) ) > 0
    % Swap h1 and h2
    h1_orig = h1;
    h2_orig = h2;
    h1( hTmp == 1 ) = h2_orig( hTmp == 1 );
    h2( hTmp == 1 ) = h1_orig( hTmp == 1 );
  end

  roots1 = quadRoots( 1, g1, h1 );
  roots2 = quadRoots( 1, g2, h2 );

  roots = [ roots1 roots2 ];
end
