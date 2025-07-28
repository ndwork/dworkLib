
function [ wavSplit, binSizes ] = makeWavSplit( sData, varargin )
  % Make the wav split to be used with the wavelet transforms
  %
  % wavSplit = makeWavSplit( sData [, minSplitSize ] )
  %
  % Written by Nicholas Dwork, Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp( 'Usage: wavSplit = makeWavSplit( sData [, minSplitSize ] )' );
    if nargout > 0, wavSplit = []; end
    return
  end

  p = inputParser;
  p.addOptional( 'minSplitSize', 16, @ispositive );
  p.parse( varargin{:} );
  minSplitSize = p.Results.minSplitSize;

  nDims = numel( sData );
  nPows = zeros( 1, nDims );

  if isscalar( minSplitSize )
    minSplitSize = minSplitSize * ones( nDims, 1 );
  end

  binSizes = zeros( 1, nDims );
  for dimIndx = 1 : nDims
    [ thisNPow, binSize ] = findNPows( sData( dimIndx ), minSplitSize( dimIndx ) );
    nPows( dimIndx ) = thisNPow;
    binSizes( dimIndx ) = binSize;
  end

  if nDims == 1
    wavSplit = zeros( 2^(nPows-1), 1 );
  else
    wavSplit = zeros( 2.^(nPows-1) );
  end
  wavSplit(1) = 1;

  sWavSplit = size( wavSplit );
  wavSplit = wavSplit( 1:min(sWavSplit), 1 : min(sWavSplit) );
end


function [ nPow, sizeDim ] = findNPows( sizeDim, minSplitSize )
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



