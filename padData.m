
function padded = padData( data, N, varargin )
  % padded = padData( data, N [, padValue ] )
  %
  % Inputs:
  % data - an array to be padded
  % N - specifies the final size of the padded data
  %
  % Optional Inputs:
  % padValue - the value of the new data elements
  %
  % Outputs:
  % padded - the padded array
  %
  % Written by Nicholas Dwork - Copyright 2016

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
    N = N * ones( numel(size(data)), 1 );
  end

  if nDimData == 1
    if isrow( data )
      padded = ones( 1, N ) * padValue;
    else
      padded = ones( N, 1 ) * padValue;
    end
    nData = numel( data );

    minY = pdFindMinIndx( N, nData );
    padded( minY : minY+nData-1 ) = data;

  elseif nDimData == 2
    padded = zeros( N(:)' );
    sData = size( data );
    minY = pdFindMinIndx( N(1), sData(1) );
    minX = pdFindMinIndx( N(2), sData(2) );
    padded( minY : minY+sData(1)-1, ...
            minX : minX+sData(2)-1 ) = data;

  elseif nDimData == 3
    padded = zeros( N(:)' );
    sData = size( data );
    minY = pdFindMinIndx( N(1), sData(1) );
    minX = pdFindMinIndx( N(2), sData(2) );
    minZ = pdFindMinIndx( N(3), sData(3) );
    padded ( minY : minY + sData(1)-1, ...
             minX : minX + sData(2)-1, ...
             minZ : minZ + sData(3)-1  ) = data;
  end

end



function minIndx = pdFindMinIndx( N, nData )
  if mod(nData,2)==0
    if mod(N,2)==0
      minIndx = (N-nData)/2 + 1;
    else
      minIndx = ceil( N/2 - nData/2 );
    end
  else
    if mod(N,2)==0
      minIndx = ceil( N/2 - nData/2 + 1 );
    else
      minIndx = (N-nData)/2 + 1;
    end
  end
end

