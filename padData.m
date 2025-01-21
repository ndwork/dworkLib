
function padded = padData( data, N, varargin )
  % padded = padData( data, N [, padValue, 'circ', trueFalse ] )
  %
  % Inputs:
  % data - an array to be padded
  % N - specifies the final size of the padded data (either a scalar, indicating
  %     the same size for all dimensions, or an array indicating size of each dimension)
  %
  % Optional Inputs:
  % padValue - the value of the new data elements (default is 0)
  %
  % Outputs:
  % padded - the padded array
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 2
    disp( 'Usage:  padded = padData( data, N [, padValue ] )' );
    return;
  end

  defaultPadValue = 0;
  p = inputParser;
  p.addOptional( 'padValue', [] );
  p.addParameter( 'circ', false );
  p.parse( varargin{:} );
  padValue = p.Results.padValue;
  circ = p.Results.circ;

  sData = size(data);
  if any( sData > N ), error('Padded size is too small.'); end

  if numel( padValue ) > 0 && circ == true
    error( 'Cannot specify a padding value and circular padding' );
  end
  if numel( padValue ) == 0
    padValue = defaultPadValue;
  end

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

    minY = findMinIndx( N, nData );
    padded( minY : minY+nData-1 ) = data;

    if circ == true
      padded( 1 : minY-1 ) = padded( nData+1 : minY+nData-1 );
      padded( minY+nData : end ) = padded( minY : N-nData );
    end

  elseif nDimData == 2
    padded = zeros( N(:)' );
    sData = size( data );
    minY = findMinIndx( N(1), sData(1) );
    minX = findMinIndx( N(2), sData(2) );
    padded( minY : minY+sData(1)-1, ...
            minX : minX+sData(2)-1 ) = data;

    if circ == true
      padded( 1 : minY-1, : ) = padded( sData(1)+1 : minY+sData(1)-1, : );
      padded( minY+sData(1) : end, : ) = padded( minY : N(1)-sData(1), : );

      padded( :, 1 : minX-1 ) = padded( :, sData(2)+1 : minX+sData(2)-1 );
      padded( :, minX+sData(2) : end ) = padded( :, minX : N(2)-sData(2) );
    end

  elseif nDimData == 3
    padded = zeros( N(:)' );
    sData = size( data );
    minY = findMinIndx( N(1), sData(1) );
    minX = findMinIndx( N(2), sData(2) );
    minZ = findMinIndx( N(3), sData(3) );
    padded ( minY : minY + sData(1)-1, ...
             minX : minX + sData(2)-1, ...
             minZ : minZ + sData(3)-1  ) = data;

    if circ == true
      padded( 1 : minY-1, :, : ) = padded( sData(1)+1 : minY+sData(1)-1, :, : );
      padded( minY+sData(1) : end, :, : ) = padded( minY : N(1)-sData(1), :, : );

      padded( :, 1 : minX-1, : ) = padded( :, sData(2)+1 : minX+sData(2)-1, : );
      padded( :, minX+sData(2) : end, : ) = padded( :, minX : N(2)-sData(2), : );

      padded( :, :, 1 : minZ-1 ) = padded( :, :, sData(3)+1 : minZ+sData(3)-1 );
      padded( :, :, minZ+sData(3) : end ) = padded( :, :, minZ : N(3)-sData(3) );
    end

  else
    padded = zeros( N(:)' );
    sData = size( data );
    minIndxs = zeros( numel(N), 1 );
    str2eval = ' padded ( ';

    if circ == true
      error( 'I still need to implement this' );
    end

    nDimsData = ndims( data );
    for dimIndx = 1 : nDimsData
      minIndxs( dimIndx ) = findMinIndx( N( dimIndx ), sData( dimIndx ) );
      str2eval = [ str2eval, ...
        ' minIndxs( ', indx2str(dimIndx,nDimsData), ...
          ' ) : minIndxs( ', indx2str(dimIndx,nDimsData), ...
        ' ) + sData( ', indx2str(dimIndx,nDimsData), ' ) - 1' ];   %#ok<AGROW>
      if dimIndx < nDimsData, str2eval = [ str2eval, ',' ]; end   %#ok<AGROW> 
    end
    str2eval = [ str2eval, ') = data;' ];
    eval( str2eval );
  end

end


function minIndx = findMinIndx( N, nData )
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

