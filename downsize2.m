
function out = downsize2( img, newSize, varargin )
  % Does nearest neighbor interpolation for image downsizing

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
