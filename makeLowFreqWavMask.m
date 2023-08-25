
function mask = makeLowFreqWavMask( sImg, varargin )
  % mask = makeLowFreqWavMask( sImg [, split ] )
  %
  % make an array of size sImg with 1s only in the lowest frequency bins of the wavelet transform
  % with recursion specified by split
  %
  % Inputs:
  % sImg - two element array specifying the size of the image
  %
  % Optional Inputs:
  % split - array specifying the number of levels of the wavelet transform.
  %   by default, split is 1 (indicating only one level).
  %   Example: [1 1 1 0] will have 3 levels.  The size of the last portion
  %   in the final level will be double that of the other portions since it
  %   wasn't split.
  %
  % Outputs:
  % mask - an array of size sImg with 1s int he lowest frequency bins and zeros everywhere else
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  defaultSplit = makeWavSplit( sImg );
  p = inputParser;
  p.addOptional( 'split', defaultSplit );
  p.parse( varargin{:} );
  split = p.Results.split;
  if numel( split ) == 0, split = defaultSplit; end

  mask = zeros( sImg );
  mask = wavMask( mask, split );
end


function out = wavMask( mask, split )

  sMask = size( mask );
  m11 = mask( 1 : sMask(1)/2, 1 : sMask(2)/2 );

  nSplit = numel(split);
  if nSplit > 1
    sSplit = size(split);
    s11 = split( 1:sSplit(1)/2, 1:sSplit(2)/2 );

    if sum( s11(:) ) > 0
      if max( mod( size( m11 ), 2 ) ) > 0
        error('wavShow: improper dimensions of image');
      end
      m11 = wavMask( m11, s11 );
    else
      m11(:) = 1;
    end
  else
    m11(:) = 1;
  end

  out = mask;
  out( 1 : sMask(1)/2, 1 : sMask(2)/2 ) = m11;
end


function wavSplit = makeWavSplit( sData, varargin )
  % Make the wav split to be used with the wavelet transforms

  p = inputParser;
  p.addOptional( 'minSplitSize', 8, @ispositive );
  p.parse( varargin{:} );
  minSplitSize = p.Results.minSplitSize;

  nDims = numel( sData );
  nPows = zeros( 1, nDims );
  for dimIndx = 1 : nDims
    nPows( dimIndx ) = findNPows( sData( dimIndx ), minSplitSize );
  end

  if nDims == 1
    wavSplit = zeros( 2^(nPows-1), 1 );
  else
    wavSplit = zeros( 2^(nPows-1) );
  end
  wavSplit(1) = 1;
end


function nPow = findNPows( sizeDim, minSplitSize )
  binPow = logBase( sizeDim, 2 );

  nPow = 0;
  for powIndx = 1 : floor( binPow )
    if mod( sizeDim, 2 ) == 0 && sizeDim/2 > minSplitSize
      nPow = nPow + 1;
      sizeDim = sizeDim / 2;
    else
      break
    end
  end

end



