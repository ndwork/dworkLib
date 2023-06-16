
function imH = showImageCube( cube, varargin )
  % imH = showImageCube( cube [, scale, 'border', border, 'borderValue', borderValue, ...
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
  %   If set to 'max' then the maximum value of the cube and the range is specified
  %     as the border value (making the border white with grayscale)
  %   If set to 'mean' then the mean value of the cube is specified as the
  %     border value.
  % scale - the amount to scale each image for display (default is 1)
  % sdevscale - passed into underlying imshowscale call
  %   (Not used if range is specified)
  % nImgsPerRow - ( default is ceil(sqrt(K)) )
  %
  % Outputs:
  % imH - a handle to the image object
  %
  % Written by Nicholas Dwork - Copyright 2018
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp([ 'Usage: imH = showImageCube( cube [, scale, ''border'', border, ''borderValue'', ', ...
      'borderValue, ''nImgsPerRow'', nImgsPerRow, ''range'', range, ''sdevScale'', sdevScale ' ]);
    return;
  end
  
  p = inputParser;
  p.addOptional( 'scale', 1, @isnumeric );
  p.addParameter( 'border', 1, @isnumeric );
  p.addParameter( 'borderValue', 'max', @(x) true);
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

  cube = squeeze( cube );

  sCube = size( cube );
  if ismatrix( cube )
    nImgs = 1;
  else
    nImgs = sCube(3);
  end
  if numel( nImgsPerRow ) == 0
    ar = 1.5;  % desired aspect ratio
    nImgsPerRow = ceil( sqrt( nImgs * sCube(1) * ar / sCube(2) ) );
  end

  if scale <= 0, error('scale must be positive'); end
  
  if scale ~= 1
    newCube = zeros( scale * sCube(1), scale * sCube(2), nImgs );
    for imgIndx = 1 : nImgs
      newCube(:,:,imgIndx) = imresize( cube(:,:,imgIndx), scale, 'nearest' );
    end
    varargin{1} = 1;
    showImageCube( newCube, varargin{:} );
    return;
  end
  
  if nImgsPerRow < 1
    error('nImgsPerRow must be positive integer');
  end
  if mod(nImgsPerRow,1) ~= 0, error('nImgsPerRow must be positive integer'); end

  nSubCols = min( nImgsPerRow, nImgs );
  nSubRows = ceil( nImgs / nImgsPerRow );

  maxCube = max( real( cube(:) ) );
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
        outImg( outImg == tmpBorderValue ) = max([ maxCube; range(:); ]);
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
      outImg( outImg == tmpBorderValue ) = max([ maxCube; range(:); ]);
    end
  end

  if numel( range ) == 0
    imH = imshowscale( outImg, scale );
  elseif isnumeric( range )
    imH = imshowscale( outImg, scale, 'range', range, 'sdevScale', sdevScale );
  else
    thresh = 0.05;
    lowScalingLevel = findValueBelowFraction( cube(:), 1-thresh );
    highScalingLevel = findValueBelowFraction( cube(:), thresh );
    scaling = [ lowScalingLevel highScalingLevel ];
  
    imH = imshowscale( outImg, scale, 'range', scaling, 'sdevScale', sdevScale );
  end

  drawnow;
end
