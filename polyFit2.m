
function c = polyFit2( x, y, z, xOrder, yOrder, varargin )
  % c = polyFit2( x, y, z, xOrder, yOrder [, 'w', w, 'cMask', cMask ] )
  %
  % Determines the matrix c that minimizes
  %    || z - \sum_{u=0,v=0}^{xOrder,yOrder} c_{u,v} x^u y^v ||_2
  %
  % Inputs:
  % x - a 1D array specifying the x coordinates
  % y - a 1D array specifying the y coordinates
  % z - a 1D array specifying p(x,y)
  % xOrder - the order of the polynomial with respect to the x values
  % yOrder - the order of the polynomial with respect to the y values
  %
  % Optional Inputs:
  % w - a 1D array specifying the weights of a weighted least squares norm
  % cMask - specify which coefficients should not be used.  It's a 2D array
  %   of size yOrder x xOrder; any value equal to 0 is not used.
  %
  % Outputs:
  % c - a 2D array of size (yOrder+1) x (xOrder+1)
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addParameter( 'w', [], @isnonnegative );
  p.addParameter( 'cMask', [], @isnumeric );
  p.parse( varargin{:} );
  w = p.Results.w;
  cMask = p.Results.cMask;

  nPts = numel( x );
  if numel( y ) ~= nPts || numel( z ) ~= nPts
    error( 'x, y, and z must have the same number of elements' );
  end

  if numel( cMask ) == 0, cMask = ones( xOrder+1, yOrder+1 ); end

  A = zeros( nPts, (xOrder+1) * (yOrder+1) );

  thisxtou = ones(nPts,1);
  for u = 0 : xOrder

    thisytov = ones(nPts,1);
    for v = 0 : yOrder

      if cMask(v+1,u+1) == 0, continue; end

      A( :, v+1 + (yOrder+1)*u ) = thisxtou .* thisytov;
      thisytov = thisytov .* y(:);

    end

    thisxtou = thisxtou .* x(:);
  end

  if numel( w ) == 0
    %c = lsqminnorm( A, z(:) );
    %c = pinv(A) * z(:);
    c = A \ z(:);
  else
    W = sparse( 1:numel(w), 1:numel(w), w(:) ); 
    %c = lsqminnorm( diag(w(:)) * A, w(:) .* z(:) );
    %c = pinv( W * A ) * ( w(:) .* z(:) );
    c = ( W * A ) \ ( w(:) .* z(:) );
  end

  c = reshape( c, [yOrder+1 xOrder+1] );
end
