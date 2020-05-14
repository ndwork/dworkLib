
function out = downsize2( img, newSize, varargin )
  % Does nearest neighbor interpolation for image downsizing

  p = inputParser;
  p.addOptional( 'op', 'notransp', @(x) true );
  p.parse( varargin{:} );
  op = p.Results.op;

  sImg = size( img );

  y = ( 0 : sImg( 1 ) - 1 ) / sImg( 1 );
  x = ( 0 : sImg( 2 ) - 1 ) / sImg( 2 );

  yNew = ( 0 : newSize( 1 ) - 1 ) / newSize( 1 );
  xNew = ( 0 : newSize( 2 ) - 1 ) / newSize( 2 );

  xMin = min( x(:) );  xMax = max( x(:) );  xNew = min( max( xNew, xMin ), xMax );
  yMin = min( y(:) );  yMax = max( y(:) );  yNew = min( max( yNew, yMin ), yMax );

  [xNew,yNew] = meshgrid( xNew, yNew );

  if strcmp( op, 'notransp' )  
    out = interp2( x, y, img, xNew, yNew, 'nearest' );

  elseif strcmp( op, 'transp' )
    out = zeros( newSize );
    yNewIndxs = round( y * newSize( 1 ) );
    xNewIndxs = round( x * newSize( 2 ) );
    out( yNewIndxs+1, xNewIndxs+1 ) = img;

  else
    error('Unrecognized op value' );
  end
end
