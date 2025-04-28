
function out = estimateDerivative( in, varargin )
  % out = estimateDerivative( in [, order] )
  %
  % estimates the first or second derivative with a central difference estimate (except at the ends of the vector)
  %
  % Written by Nicholas Dwork - Copyright 2025
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addOptional( 'order', 1, @(x) x==1 || x==2 );
  p.addParameter( 'd', 1, @ispositive );
  p.parse( varargin{:} );
  order = p.Results.order;
  d = p.Results.d;

  out = zeros( size( in ) );

  if order == 1

    if numel( in ) < 2, error( 'in does not have enough elements' ); end

    out(1) = in(2) - in(1);
    out(2:end-1) = ( in(3:end) - in(1:end-2) ) / (2*d);
    out(end) = in(end) - in(end-1);

  elseif order == 2

    if numel( in ) < 3, error( 'in does not have enough elements' ); end

    out(1) = ( in(3) - 2*in(2) + in(1) ) / (d*d);
    out(2:end-1) = ( in(3:end) - 2*in(2:end-1) + in(1:end-2) ) / (d*d);
    out(end) = ( in(end) - 2*in(end-1) + in(end-2) ) / (d*d);

  else
    error( 'Unrecognized value of order' );
  end

end
