
function recon = mri_reconStructuredSparseSENSE( kData, sMaps, lambda, varargin )
  % recon = mri_reconStructuredSparseSENSE( kData, sMaps, lambda, [, 'img0', img0, 'nIter', nIter, ...
  %   'noiseCov', noiseCov, 'transformType', transformType, 'wavSplit', wavSplit ] )
  %
  % Inputs:
  % kData - an array of size Ny x Nx x nCoils
  %
  % Outputs:
  % recon - a 2D complex array that is the image
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  sImg = [ size( kData, 1 ) size( kData, 2 ) ];

  wavSplit = [];
  transformType = [];
  for vIndx = 1 : numel( varargin )-1
    if strcmp( varargin{vIndx}, 'wavSplit' )
      wavSplit = varargin{ vIndx + 1 };
    end
    if strcmp( varargin{vIndx}, 'transformType' )
      transformType = varargin{ vIndx + 1 };
    end
  end

  if numel( wavSplit ) == 0, wavSplit = makeWavSplit( sImg ); end
  if numel( transformType ) == 0, transformType = 'wavelet'; end

  sImg = size( kData );
  if strcmp( transformType, 'curvelet' )
    acr = makeCurvAutoCalRegion( sImg );

  elseif strcmp( transformType, 'wavelet' )
    acr = makeWavAutoCalRegion( sImg, wavSplit );

  elseif strcmp( transformType, 'wavCurv' )
    acrCurv = makeCurvAutoCalRegion( sImg );
    acrWav = makeWavAutoCalRegion( sImg, wavSplit );
    acr = max( acrCurv, acrWav );

  else
    error( 'Unrecognized transform type' );
  end

  acrDataL = bsxfun( @times, acr, kData );
  beta = kData - acrDataL;
  beta( kData == 0 ) = 0;

  coilReconsL = mri_reconIFFT( acrDataL, 'multiSlice', true );
  reconL = mri_reconRoemer( coilReconsL );

  reconH = mri_reconSparseSENSE( beta, sMaps, lambda );

  recon = reconH + reconL;
end



function acr = makeCurvAutoCalRegion( sImg, varargin )
  
  p = inputParser;
  p.addOptional( 'beta', 4, @ispositive );
  p.parse( varargin{:} );
  beta = p.Results.beta;

  x = zeros( sImg );
  [~,lowPass] = fdct_wrapping( x, false, 1 );

  sLowPass = size( lowPass );
  kx = kaiser( sLowPass(2), beta );
  ky = kaiser( sLowPass(1), beta );
  [ kxs, kys ] = meshgrid( kx, ky );

  acr = zeros( sImg );
  acr(1:sLowPass(1),1:sLowPass(2)) = kxs .* kys;

  acr = circshift( acr, [ -floor(sLowPass(1) * 0.5) -floor(sLowPass(2) * 0.5) ] );
  acr = fftshift( acr );
end


function acr = makeWavAutoCalRegion( sImg, wavSplit, varargin )

  p = inputParser;
  p.addOptional( 'beta', 4, @ispositive );
  p.parse( varargin{:} );
  beta = p.Results.beta;

  wavMask = makeLowFreqWavMask( sImg, wavSplit );

  lastX = find( wavMask(1,:), 1, 'last' );
  lastY = find( wavMask(:,1), 1, 'last' );

  kx = kaiser( lastX, beta );
  ky = kaiser( lastY, beta );
  [ kxs, kys ] = meshgrid( kx, ky );

  acr = zeros( size( wavMask ) );
  acr(1:lastY,1:lastX) = kxs .* kys;
  %acr(1:lastY,1:lastX) = 1;

  acr = circshift( acr, [ -floor(lastY * 0.5) -floor(lastX * 0.5) ] );
  acr = fftshift( acr );
end


function wavSplit = makeWavSplit( sImg )

  minSplit = 8;

  yTmp = sImg(1);
  ySplitSize = 1;
  while ( mod( yTmp, 1 ) == 0 )
    ySplitSize = ySplitSize * 2;
    yTmp = yTmp / 2;
  end
  ySplitSize = ySplitSize / 4;
  ySplitSize = min( ySplitSize, minSplit );

  xTmp = sImg(2);
  xSplitSize = 1;
  while ( mod( xTmp, 1 ) == 0 )
    xSplitSize = xSplitSize * 2;
    xTmp = xTmp / 2;
  end
  xSplitSize = xSplitSize / 4;
  xSplitSize = min( xSplitSize, minSplit );

  wavSplit = zeros( [ ySplitSize xSplitSize ] );
  wavSplit(1) = 1;
end


