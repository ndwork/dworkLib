
function showImageCube( cube, varargin )
  % showImageCube( cube [, scale, 'nImgsPerRow', nImgsPerRow ] )
  % Given an MxNxK array, shows each slice as a separate image on a single
  % figure.
  %
  % Inputs:
  % cube is an MxNxK array
  %
  % Optional Inputs:
  % scale - the amount to scale each image for display (default is 1)
  % nImgsPerRow - (default is 5)
  %
  % Written by Nicholas Dwork - Copyright 2018
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addOptional( 'scale', 1, @isnumeric );
  p.addParameter( 'nImgsPerRow', 5, @isnumeric );
  p.parse( varargin{:} );
  scale = p.Results.scale;
  nImgsPerRow = p.Results.nImgsPerRow;

  if scale <= 0, error('scale must be positive'); end;
  if nImgsPerRow <= 1, error('nImgsPerRow must be positive integer'); end;
  if mod(nImgsPerRow,1) ~= 0, error('nImgsPerRow must be positive integer'); end;

  sCube = size( cube );
  nImgs = sCube(3);
  nSubCols = min( nImgsPerRow, nImgs );
  nSubRows = ceil( nImgs / nImgsPerRow );
  
  outImg = inplaceImg( cube(:,:,1), nSubRows, nSubCols, 1 );
  for i=2:nImgs
    outImg = inplaceImg( cube(:,:,i), nSubRows, nSubCols, i, outImg );
  end

  figure; imshowscale( outImg, scale );  drawnow;
end
