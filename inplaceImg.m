
function outImg = inplaceImg( subImg, nSubRows, nSubCols, subIndx, varargin )
  % outImg = inplaceImg( subImg, nSubRows, nSubCols, subIndx [, inImg, ...
  %   'order', order, 'border', border, 'borderValue', borderValue ] )
  %
  % Inputs:
  % subImg - the image to be put into the out image
  % nSubRows - the number of image rows in outImg
  % nSubCols - the number of image cols in outImg
  % subIndx - the (one based) index of the current image to input
  %
  % Optional Inputs:
  % inImg - the input image to be altered (if not provided, creates an
  %   image of all zeros)
  % order - either 'rowMajor' (default) or 'colMajor'
  % border - specifies a border to place between images (default is 0)
  % borderValue - specifies the value to place in between images.
  %   (Only used when inImg is not specified.)
  %
  % Outputs:
  % outImg - the complete image
  %
  % Ex:
  % im = phantom();
  % completeImg = inplaceImg( im, 2, 3, 1 );
  % im = flipud( im );
  % completeImg = inplaceImg( im, 2, 3, 2, completeImg );
  % im = fliplr( im );
  % completeImg = inplaceImg( im, 2, 3, 3, completeImg );
  %
  % Written by Nicholas Dwork - Copyright 2017
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  defaultBorder = 0;
  defaultBorderValue = 0;
  p = inputParser;
  p.addOptional( 'inImg', [] );
  p.addParameter( 'order', 'rowMajor' );
  p.addParameter( 'border', defaultBorder, @isnumeric );
  p.addParameter( 'borderValue', defaultBorderValue, @isnumeric );
  p.parse( varargin{:} );
  border = p.Results.border;
  borderValue = p.Results.borderValue;
  inImg = p.Results.inImg;
  order = p.Results.order;

  nRows = size( subImg, 1 );
  nCols = size( subImg, 2 );
  
  if numel( inImg ) == 0
    outImg = double( borderValue ) * ...
      ones( nSubRows*nRows + border*(nSubRows-1), ...
            nSubCols*nCols + border*(nSubCols-1), size(subImg,3) );
  else
    outImg = inImg;
  end
  sOutImg = size( outImg );

  if strcmp( order, 'colMajor' )
    subRowIndx = mod( subIndx-1, nSubRows ) + 1;
    subColIndx = ceil( subIndx / nSubRows );
  else
    subColIndx = mod( subIndx-1, nSubCols ) + 1;
    subRowIndx = ceil( subIndx / nSubCols );
  end

  colIndxL = (subColIndx-1) * nCols + (subColIndx-1)*border + 1;
  colIndxH = min( subColIndx * nCols + (subColIndx-1)*border, sOutImg(2) );

  rowIndxL = (subRowIndx-1) * nRows + (subRowIndx-1)*border + 1;
  rowIndxH = min( subRowIndx * nRows + (subRowIndx-1)*border, sOutImg(1) );

  dCol = colIndxH - colIndxL + 1;
  dRow = rowIndxH - rowIndxL + 1;
  if dCol > 0 && dRow > 0
    outImg( rowIndxL:rowIndxH, colIndxL:colIndxH, : ) = subImg(1:dRow,1:dCol,:);
  end
end

