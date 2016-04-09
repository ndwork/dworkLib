
function padded = padData( data, N, varargin )
  % padded = padData( data, N [, padValue ] )
  % data - an array to be padded
  % N - specifies the final size of the padded data

  defaultPadValue = 0;
  p = inputParser;
  p.addOptional( 'padValue', defaultPadValue );
  p.parse( varargin{:} );
  padValue = p.Results.padValue;

  sData = size(data);
  tooSmall = max( sData > N );
  if tooSmall>0, error('Padded size is too small.'); end;

  nDimData = ndims( data );
  if nDimData==2 && min( size(data) ) == 1
    nDimData = 1;
  end

  if numel(N) ~= nDimData
    error('Incorrect size of padded data specified');
  end

  if nDimData == 1
    if isrow( data )
      padded = ones( 1, N ) * padValue;
    else
      padded = ones( N, 1 ) * padValue;
    end
    nData = numel( data );
    minY = ceil( N/2 - nData/2 + 1 );
    padded( minY : minY+nData-1 ) = data;

  elseif nDimData == 2
    padded = zeros( N );
    sData = size( data );
    minY = ceil( N(1)/2 - sData(1)/2 + 1 );
    minX = ceil( N(2)/2 - sData(2)/2 + 1 );
    padded( minY : minY+sData(1)-1, ...
            minX : minX+sData(2)-1 ) = data;

  elseif nDimData == 3
    padded = zeros( N );
    sData = size( data );
    minY = ceil( N(1)/2 - sData(1)/2 + 1 );
    minX = ceil( N(2)/2 - sData(2)/2 + 1 );
    minZ = ceil( N(3)/2 - sData(3)/2 + 1 );
    padded ( minY : minY + sData(1)-1, ...
             minX : minX + sData(2)-1, ...
             minZ : minZ + sData(3)-1  ) = data;
  end


end
