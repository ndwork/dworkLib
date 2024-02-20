

function wavSplit = makeWavSplit( sData, varargin )
  % Make the wav split to be used with the wavelet transforms

  p = inputParser;
  p.addOptional( 'minSplitSize', 16, @ispositive );
  p.parse( varargin{:} );
  minSplitSize = p.Results.minSplitSize;

  nDims = numel( sData );
  nPows = zeros( 1, nDims );

  if numel( minSplitSize ) == 1
    minSplitSize = minSplitSize * ones( nDims, 1 );
  end

  for dimIndx = 1 : nDims
    nPows( dimIndx ) = findNPows( sData( dimIndx ), minSplitSize( dimIndx ) );
  end

  if nDims == 1
    wavSplit = zeros( 2^(nPows-1), 1 );
  else
    wavSplit = zeros( 2.^(nPows-1) );
  end
  wavSplit(1) = 1;
end


function nPow = findNPows( sizeDim, minSplitSize )
  binPow = logBase( sizeDim, 2 );

  nPow = 0;
  for powIndx = 1 : floor( binPow )
    if mod( sizeDim, 2 ) == 0 && sizeDim/2 >= minSplitSize
      nPow = nPow + 1;
      sizeDim = sizeDim / 2;
    else
      break
    end
  end

end



