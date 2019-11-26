
function p = evaluatePoly2( c, x, y )
  % p = evaluatePoly2( c, x, y )
  %
  % This function evaluates the 2D polynomial defined by coefficients in c
  % at points (x,y).
  %
  % Inputs:
  % x - 1D array of values
  % y - 1D array of values
  %
  % Outputs:
  % p = a 2 dimensional array of polynomial coefficients.
  %   p(x,y) = \sum_{u=0,v=0}^{xOrder,yOrder} c_{u,v} x^u y^v
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  nPts = numel( x );
  if numel(y) ~= nPts, error( 'x and y must have the same number of elements' ); end

  xOrder = size( c, 2 ) - 1;
  yOrder = size( c, 1 ) - 1;

  p = zeros( nPts, 1 );

  xtou = ones( nPts, 1 );
  for u = 0 : xOrder

    ytov = ones( nPts, 1 );
    for v = 0 : yOrder
      p = p + c(v+1,u+1) * xtou .* ytov;
      ytov = ytov .* y(:);
    end

    xtou = xtou .* x(:);
  end

  p = reshape( p, size(x) );
end
