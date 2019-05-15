
function showImageCube( cube, varargin )
  % showImageCube( cube [, scale, 'border', border, 'borderValue', borderValue, ...
  %   'nImgsPerRow', nImgsPerRow, 'range', range, 'sdevScale', sdevScale ] )
  %
  % Given an MxNxK array, shows each slice as a separate image on a single
  % figure.
  %
  % Inputs:
  % cube is an MxNxK array
  %
  % Optional Inputs:
  % border - specifies a border to place between images (default is 0)
  % borderValue - specifies the value to place in between images
  %   If set to 'max' then the maximum value of the cube is specified as the
  %     border value (making the border white with grayscale)
  %   If set to 'mean' then the mean value of the cube is specified as the
  %     border value.
  % scale - the amount to scale each image for display (default is 1)
  % sdevscale - passed into underlying imshowscale call
  %   (Not used if range is specified)
  % nImgsPerRow - ( default is ceil(sqrt(K)) )
  %
  % Written by Nicholas Dwork - Copyright 2018
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addOptional( 'scale', 1, @isnumeric );
  p.addParameter( 'border', 0, @isnumeric );
  p.addParameter( 'borderValue', 0, @(x) true);
  p.addParameter( 'nImgsPerRow', [], @isnumeric );
  p.addParameter( 'range', [] );
  p.addParameter( 'sdevScale', [], @isnumeric );
  p.parse( varargin{:} );
  scale = p.Results.scale;
  border = p.Results.border;
  borderValue = p.Results.borderValue;
  nImgsPerRow = p.Results.nImgsPerRow;
  range = p.Results.range;
  sdevScale = p.Results.sdevScale;

  sCube = size( cube );
  if ismatrix( cube )
    nImgs = 1;
  else
    nImgs = sCube(3);
  end
  if numel( nImgsPerRow ) == 0, nImgsPerRow = ceil( sqrt( nImgs ) ); end

  if scale <= 0, error('scale must be positive'); end
  if nImgsPerRow < 1
    error('nImgsPerRow must be positive integer');
  end
  if mod(nImgsPerRow,1) ~= 0, error('nImgsPerRow must be positive integer'); end

  nSubCols = min( nImgsPerRow, nImgs );
  nSubRows = ceil( nImgs / nImgsPerRow );

  maxCube = max( cube(:) );
  tmpBorderValue = maxCube + 1;
  outImg = inplaceImg( cube(:,:,1), nSubRows, nSubCols, 1, ...
    'border', border, 'borderValue', tmpBorderValue );
  for i=2:nImgs
    outImg = inplaceImg( cube(:,:,i), nSubRows, nSubCols, i, outImg, ...
    'border', border, 'borderValue', tmpBorderValue );
  end
  if border > 0
    if isnumeric( borderValue )
      outImg( outImg == tmpBorderValue ) = borderValue;
    else
      if strcmp( borderValue, 'max' )
        outImg( outImg == tmpBorderValue ) = maxCube;
      elseif strcmp( borderValue, 'mean' )
        outImg( outImg == tmpBorderValue ) = mean( cube(:) );
      else
        error( 'Unrecognized border value term' );
      end
    end
  else
    if isnumeric( borderValue )
      outImg( outImg == tmpBorderValue ) = borderValue;
    else
      outImg( outImg == tmpBorderValue ) = maxCube;
    end
  end

  imshowscale( outImg, scale, 'range', range, 'sdevScale', sdevScale );
  drawnow;
end
