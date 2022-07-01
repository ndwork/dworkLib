
function out = circularCrop( data, newSize )
  % out = circularCrop( data, newSize )
  %
  % Inputs:
  % data - array
  % newSize - 1D array with number of elements equal to the number of dimensions of data
  %           or a scalar (which assumes same cropped size in all dimensions)
  %
  % Outputs:
  % out - array of size newSize
  %
  % Written by Nicholas Dwork, Copyright 2022
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It 
  % is offered without any warranty expressed or implied, including the 
  % implied warranties of merchantability or fitness for a particular 
  % purpose.

  sData = size( data );

  if numel( sData ) == 1
    nDimsData = 1;
  elseif ( numel( sData ) == 2 )  &&  ( sData(2) == 1 )
      nDimsData = 1;
      sData = sData(1);
  else
    nDimsData = ndims( data );
  end

  if ( nDimsData ~= numel( newSize) ) && ( numel( newSize ) ~= 1 )
    error( 'Must provide a newSize for each dimension' );
  end
  if nDimsData > 1 && numel( newSize ) == 1
    newSize = newSize * ones( ndims(data), 1 );
  end
  if max( newSize(:) > sData(:) ) == 1
    error( 'New size must be <= the size of data' );
  end

  croppedCoords = cell( nDimsData, 1 );
  for dim = 1 : nDimsData
    croppedCoords{ dim } = getCroppedCoords( sData(dim), newSize(dim) );
  end

  if numel( newSize ) ~= 2, error( 'Only two dimensions implemented right now' ); end

  tmp = zeros( newSize(1), sData(2) );
  tmpCoords = size2imgCoordinates( newSize(1) );
  for i = 1 : numel( tmpCoords )
    thisCoord = tmpCoords(i);
    tmp(i,:) = sum( data( croppedCoords{1} == thisCoord, : ), 1 );
  end

  out = zeros( newSize );
  outCoords = size2imgCoordinates( newSize(2) );
  for i = 1 : numel( outCoords )
    thisCoord = outCoords(i);
    out(:,i) = sum( tmp( :, croppedCoords{2} == thisCoord ), 2 );
  end

end



function outCoords = getCroppedCoords( dataSize, newSize )

  if newSize > dataSize
    error( 'newSize must be <= dataSize' );

  elseif newSize == dataSize
    outCoords = size2imgCoordinates( newSize );

  else
    newCoords = size2imgCoordinates( newSize );
    maxNewCoord = max( newCoords );
    minNewCoord = min( newCoords );
  
    outCoords = size2imgCoordinates( dataSize );
    while true
      if max( outCoords ) <= maxNewCoord, break; end
  
      outCoords( outCoords > maxNewCoord ) = outCoords( outCoords > maxNewCoord ) - newSize;
    end
  
    while true
      if min( outCoords ) >= minNewCoord, break; end
  
      outCoords( outCoords < minNewCoord ) = outCoords( outCoords < minNewCoord ) + newSize;
    end

  end
end


