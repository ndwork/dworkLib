
function out = downsize2( img, newSize, varargin )
  % Does nearest neighbor interpolation for image downsizing
  %
  % out = downsize2( img, newSize [, op ] )
  %
  % Inputs:
  % img - 2D array
  % newSize - two element array specifying new size of image
  %
  % Optional Inputs:
  % op - either 'transp' or 'notransp'
  %
  % Outputs:
  % out - the downsized image
  %
  % Written by Nicholas Dwork - Copyright 2020
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addOptional( 'op', 'notransp', @(x) true );
  p.parse( varargin{:} );
  op = p.Results.op;

  sImg = size( img );

  y = ( ( 0 : sImg( 1 ) - 1 ) + 0.5 ) / sImg( 1 );
  x = ( ( 0 : sImg( 2 ) - 1 ) + 0.5 ) / sImg( 2 );

  yNew = ( ( 0 : newSize( 1 ) - 1 ) + 0.5 ) / newSize( 1 );
  xNew = ( ( 0 : newSize( 2 ) - 1 ) + 0.5 ) / newSize( 2 );


  if strcmp( op, 'notransp' )  
    xNewIndxs = interp1( x, 1:numel(x), xNew, 'nearest' );
    yNewIndxs = interp1( y, 1:numel(y), yNew, 'nearest' );
    out = img( yNewIndxs, xNewIndxs );

  elseif strcmp( op, 'transp' )
    xIndxs = interp1( xNew, 1:numel(xNew), x, 'nearest' );
    yIndxs = interp1( yNew, 1:numel(yNew), y, 'nearest' );
    out = zeros( newSize );
    out( yIndxs, xIndxs ) = img;

  else
    error('Unrecognized op value' );
  end
end
