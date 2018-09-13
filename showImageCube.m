
function showImageCube( cube, varargin )
  showImageCube( cube [, scale] )
  % Given an MxNxK array, shows each slice as a separate image on a single
  % figure.
  %
  % Inputs:
  % cube is an MxNxK array
  % scale (optional) is the amount to scale each image for display
  %
  % Written by Nicholas Dwork - Copyright 2018
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addOptional( 'scale', 1, @isnumeric );
  p.parse( varargin{:} );
  scale = p.Results.scale;

  sCube = size( cube );
  nImgs = sCube(3);
  nSubCols = min( 4, nImgs );
  nSubRows = ceil( nImgs / 4 );
  
  outImg = inplaceImg( cube(:,:,1), nSubRows, nSubCols, 1 );
  for i=2:nImgs
    outImg = inplaceImg( cube(:,:,1), nSubRows, nSubCols, i, outImg );
  end

  figure; imshowscale( outImg, scale );
end
