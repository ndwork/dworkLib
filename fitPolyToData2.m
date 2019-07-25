
function c = fitPolyToData2( xOrder, yOrder, x, y, z )
  % c = fitPolyToData2( xOrder, yOrder, x, y, z ) or
  % c = fitPolyToData2( xOrder, yOrder, z )
  %
  % This function finds a polynomial p so that || z - p(x,y) ||_2 is minimized
  %
  % Inputs:
  % xOrder - the order of the polynomial in x
  % yOrder - the order of the polynomial in y
  % y - (optional) domain values.  If not supplied, x = 1, 2, ..., size(z,1)
  % x - (optional) domain values.  If not supplied, x = 1, 2, ..., size(z,2)
  % z - range values
  %
  % Outputs:
  % c = a 2 dimensional array of polynomial coefficients.
  %   p(x,y) = \sum_{u=0,v=0}^{xOrder,yOrder} c_{u,v} x^u y^v
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  Nx = numel( x );
  Ny = numel( y );

  if Nx ~= numel( z ), error( 'Number of elements of x not compatible with z' ); end
  if Ny ~= numel( z ), error( 'Number of elements of y not compatible with z' ); end

  A = ones( Nx, (xOrder+1) * (yOrder+1) );
  colIndx = 1;
  for u = 0 : xOrder
    xTou = x(:).^u;

    for v = 0 : yOrder
      A(:,colIndx) = xTou .* y(:).^v;
      colIndx = colIndx + 1;
    end
  end

  c = A \ z(:);
  c = reshape( c, [ yOrder+1 xOrder+1 ] );
end

