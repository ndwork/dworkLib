
function [recon,lambda] = mri_reconStructuredSparseSENSE( kData, sMaps, lambda, varargin )
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

  p = inputParser;
  p.addParameter( 'img0', [], @isnumeric );
  p.addParameter( 'noiseCov', [], @isnumeric );
  p.addParameter( 'nReweightIter', [], @ispositive );
  p.addParameter( 'optAlg', [], @(x) true );
  p.addParameter( 'reweightEpsilon', [], @ispositive );
  p.addParameter( 't', [], @ispositive );
  p.addParameter( 'transformType', 'wavelet', @(x) true );
  p.addParameter( 'waveletType', 'Daubechies-4', @(x) true );
  p.addParameter( 'wavSplit', [], @isnumeric );
  p.parse( varargin{:} );
  img0 = p.Results.img0;
  noiseCov = p.Results.noiseCov;
  nReweightIter = p.Results.nReweightIter;
  optAlg = p.Results.optAlg;
  reweightEpsilon = p.Results.reweightEpsilon;
  t = p.Results.t;
  transformType = p.Results.transformType;
  waveletType = p.Results.waveletType;
  wavSplit = p.Results.wavSplit;

  if numel( transformType ) == 0, transformType = 'wavelet'; end

  sImg = [ size( kData, 1 ) size( kData, 2 ) ];
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

  if numel( img0 ) == 0
    coilRecons = fftshift2( ifft2( ifftshift2( kData ) ) );
    img0 = mri_reconRoemer( coilRecons, 'sMaps', sMaps );
  end

  fftImg0 = fftshift2( fft2( ifftshift2( img0 ) ) );
  fftImg0_H = fftImg0 .* ( 1 - acr );
  img0_H = fftshift2( ifft2( ifftshift2( fftImg0_H ) ) );

  acrK = bsxfun( @times, kData, acr );
  reconLs = mri_reconIFFT( acrK, 'multiSlice', true );

  beta = kData - acrK;
  beta( kData == 0 ) = 0;
  [reconH,lambda] = mri_reconSparseSENSE( beta, sMaps, lambda, 'img0', img0_H, 'noiseCov', noiseCov, ...
    'nReweightIter', nReweightIter, 'optAlg', optAlg, 'reweightEpsilon', reweightEpsilon, 't', t, ...
    'transformType', transformType, 'waveletType', waveletType, 'wavSplit', wavSplit );
  reconHs = bsxfun( @times, sMaps, reconH );

  coilRecons = reconLs + reconHs;
  recon = mri_reconRoemer( coilRecons );
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

  acr = circshift( acr, [ -floor(lastY * 0.5) -floor(lastX * 0.5) ] );
  acr = fftshift( acr );
end

