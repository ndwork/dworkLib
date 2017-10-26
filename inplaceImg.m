
function outImg = inplaceImg( subImg, nSubRows, nSubCols, subIndx, varargin )
  % outImg = inplaceImg( subImg, nSubRows, nSubCols, subIndx, varargin )
  %
  % Inputs:
  % subImg - the image to be put into the out image
  % nSubRows - the number of image rows in outImg
  % nSubCols - the number of image cols in outImg
  % subIndx - the index of the current image to input
  % inImg (optional) - the input image to be altered (if not provided, creates an
  %   image of all zeros)
  % order (optional) - either 'rowMajor' (default) or 'colMajor'
  %
  % Written by Nicholas Dwork - Copyright 2017
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addOptional( 'inImg', [] );
  p.addParameter( 'order', 'rowMajor' );
  p.parse( varargin{:} );
  inImg = p.Results.inImg;
  order = p.Results.order;

  nRows = size( subImg, 1 );
  nCols = size( subImg, 2 );

  if numel( inImg ) == 0
    outImg = zeros( nSubRows*nRows, nSubCols*nCols, size(subImg,3) );
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

  colIndxL = (subColIndx-1) * nCols + 1;
  colIndxH = min( subColIndx * nCols, sOutImg(2) );

  rowIndxL = (subRowIndx-1) * nRows + 1;
  rowIndxH = min( subRowIndx * nRows, sOutImg(1) );

  outImg( rowIndxL:rowIndxH, colIndxL:colIndxH, : ) = subImg;
end

