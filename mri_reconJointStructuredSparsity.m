
function [recon,sMaps,oValues,lambda] = mri_reconJointStructuredSparsity( kData, varargin )
  % recon = mri_reconJointStructuredSparsity( samples [, 'lambda', lambda, 'nIter', nIter, 
  %   'nReweightIter', nReweightIter, 'polyOrder', polyOrder, 'printEvery', printEvery, ...
  %   'transformType', transformType, 'verbose', verbose, 'waveletType', waveletType, ...
  %   'w', w, 'wavSplit', wavSpit ] )
  %
  % This routine uses proximal gradient methods to minimize
  %   0.5 * || A y - b ||_2^2 + lambda || y ||_{w,1}
  %   where A is sampleMask * Fourier Transform * adjoint( Psi ).
  %   Here, Psi is either a wavelet or curvelet transform.  Both are orthogonal.
  %   The reconstruction returned is adjoint( Psi ) * y.
  %
  % Note that || y ||_{w,1} = w1 |y1| + w2 |y2| + ... + wN |yN|
  %
  % Iterative reweighting is done according to "Enhancing sparsity by reweighted L1
  %   minimization" by Candes, Watkin, and Boyd with the nReweightIter optional parameter.
  %
  % Inputs:
  % samples - a 2D array that is zero wherever a sample wasn't acquired
  % lambda - regularization parameter
  %
  % Optional Inputs:
  % alg - choose the algorithm to perform the minimization.  Default is 'FISTA_wLS'.
  %   Options are FISTA_wLS.
  % nIter - the number of iterations that FISTA will perform (default is 100)
  % nReweightIter - number of reweighting iterations
  % printEvery - FISTA prints a verbose statement every printEvery iterations
  % transformType - Choose the sparifying transformation.  Default is WavCurv (which uses
  %   the wavelet transform and the curvelet transform as a redundant dictionary).
  %   Options are:  'curvlet', 'wavelet', or 'wavCurv'.
  %     curvelet - Discrete curvelet transform (requires CurveLab; see curvelet.org)
  %     wavelet - Wavelet transform (type specified by waveletType parameter)
  %     wavCurv - Redundant dictionary of wavelet and curvelet (default)
  % verbose - if true, prints informative statements
  % waveletType - either 'Daubechies-4' (default) or 'Haar'
  %
  % Written by Nicholas Dwork - Copyright 2023
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  sData = size( kData );
  sImg = sData(1:2);
  nImg = prod( sImg );
  nCoils = sData( 3 );
  wavSplit = makeWavSplit( sImg );

  p = inputParser;
  p.addParameter( 'lambda', [], @isnumeric );
  p.addParameter( 'alg', 'fista_wLS', @(x) true );
  p.addParameter( 'checkAdjoints', false, @islogical );
  p.addParameter( 'nIter', [], @ispositive );
  p.addParameter( 'noiseCov', [], @isnumeric );
  p.addParameter( 'nReweightIter', [], @ispositive );
  p.addParameter( 'polyOrder', [], @ispositive );
  p.addParameter( 'printEvery', 1, @ispositive );
  p.addParameter( 'stepSize', [], @ispositive );
  p.addParameter( 'transformType', 'wavelet', @(x) true );
  p.addParameter( 'verbose', false, @(x) isnumeric(x) || islogical(x) );
  p.addParameter( 'w', 1, @isnumeric );
  p.addParameter( 'waveletType', [], @(x) true );
  p.addParameter( 'wavSplit', wavSplit, @isnumeric );
  p.parse( varargin{:} );
  lambda = p.Results.lambda;
  alg = p.Results.alg;
  checkAdjoints = p.Results.checkAdjoints;
  nIter = p.Results.nIter;
  noiseCov = p.Results.noiseCov;
  nReweightIter = p.Results.nReweightIter;
  polyOrder = p.Results.polyOrder;
  printEvery = p.Results.printEvery;
  stepSize = p.Results.stepSize;
  transformType = p.Results.transformType;
  verbose = p.Results.verbose;
  w = p.Results.w;
  waveletType = p.Results.waveletType;
  wavSplit = p.Results.wavSplit;

  if numel( nReweightIter ) == 0 || nReweightIter == 0
    if numel( lambda ) == 0 || lambda == 0
      nReweightIter = 10;
    else
      nReweightIter = 1;
    end
  end
  reweightEpsilon = 0.1;

  if numel( transformType ) == 0, transformType = 'wavCurv'; end
  if numel( waveletType ) == 0, waveletType = 'Daubechies-4'; end

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

  if numel( noiseCov ) > 0
    invNoiseCov = inv( noiseCov );
    [~,s,~] = svd( invNoiseCov, 'econ' );
    invNoiseCov = invNoiseCov ./ s(1);
    L = chol( invNoiseCov, 'lower' );
  end
  function out = applyL( in, type )
    if numel( noiseCov ) == 0, out = in; return; end

    % Assumes last dimension is coil dimension
    if nargin < 2, type = 'notransp'; end

    sIn = size( in );
    reshaped = reshape( in, [ prod( sIn(1:end-1) ) sIn(end) ] );
    if strcmp( type, 'notransp' )
      out = transpose( L * transpose( reshaped ) );
    else
      out = transpose( L' * transpose( reshaped ) );
    end
    out = reshape( out, sIn );
  end

  M = ( abs( kData ) ~= 0 );
  acrKDataL = bsxfun( @times, kData, acr );
  beta = kData - acrKDataL;
  beta = beta( M == 1 );
  LStarBeta = applyL( beta, 'transp' );
  LStarBeta = LStarBeta(:);

  coilReconsL = mri_reconIFFT( acrKDataL );
  reconL = mri_reconRoemer( coilReconsL );

  function out = F( x )
    out = fftshift2( fft2( ifftshift2( x ) ) );
  end

  function out = FH( y )
    out = fftshift2( fft2h( ifftshift2( y ) ) );
  end

  coilRecons = mri_reconIFFT( kData );
  intensityMask = mri_makeIntensityMask( reshape( kData, sImg(1), sImg(2), 1, [] ), ...
    'noiseCoords', [ 1 15 1 15 ] );
  ssqRecon = mri_reconSSQ( kData );
  sMaps = bsxfun( @rdivide, coilRecons, ssqRecon );
  sMaps = bsxfun( @times, sMaps, intensityMask );

  function out = applyS( in, type )
    if nargin < 2 || strcmp( type, 'notransp' )
      out = bsxfun( @times, sMaps, in );
    else
      out = sum( conj(sMaps) .* in, 3 );
    end
  end

  function out = findCellSizes( in )
    if iscell( in )
      out = cell( size(in) );
      for indx = 1 : numel( in )
        theseCellSizes = findCellSizes( in{ indx } );
        out{ indx } = theseCellSizes;
      end
    else
      out = size( in );
    end
  end

  function nTotal = sumCurvCellSizes( cellSizes )
    nTotal = 0;
    if iscell( cellSizes )
      for indx = 1 : numel( cellSizes )
        thisCell = cellSizes{ indx };  
        nTotal = nTotal + sumCurvCellSizes( thisCell );
      end
    else
      nTotal = nTotal + prod( cellSizes );
    end
  end

  function out = curvCell2Vec( cx )
    if iscell( cx )
      out = cell( numel( cx ), 1 );
      for indx = 1 : numel( cx )
        out{ indx } = curvCell2Vec( cx{ indx } );
      end
      out = cell2mat( out );
    else
      out = cx(:);
    end
  end

  function out = vec2CurvCell( v, curvCellSizes )
    if iscell( curvCellSizes )
      out = cell( size( curvCellSizes ) );
      thisIndx = 1;
      for indx = 1 : numel( curvCellSizes )
        nSubVec = sumCurvCellSizes( curvCellSizes{ indx } );
        subVec = v( thisIndx : thisIndx + nSubVec - 1 );
        out{ indx } = vec2CurvCell( subVec, curvCellSizes{ indx } );
        thisIndx = thisIndx + nSubVec;
      end
    else
      out = reshape( v, curvCellSizes );
    end
  end

  function out = curvelet( x, type )
    if nargin < 2 || strcmp( type, 'notransp' )
      curvCells = fdct_wrapping( x, false );
      out = curvCell2Vec( curvCells );
    else
      curvCells = vec2CurvCell( x, curvCellSizes );
      out = ifdct_wrapping( curvCells, false );
    end
  end

  if strcmp( waveletType, 'Daubechies-4' )
    wavTrans = @(x) wtDaubechies2( x, wavSplit );
    wavTransH = @(y) iwtDaubechies2( y, wavSplit );
  elseif strcmp( waveletType, 'Haar' )
    wavTrans = @(x) wtHaar2( x, wavSplit );
    wavTransH = @(y) iwtHaar2( y, wavSplit );
  else
    error( 'Unrecognized wavelet type' );
  end

  function out = wavCurv( x, type )
    if nargin < 2 || strcmp( type, 'notransp' )
      x = reshape( x, sImg );
      wx = wavTrans( x );
      cx = curvelet( x );
      out = [ wx(:); cx(:); ];
    elseif strcmp( type, 'transp' )
      x1 = reshape( x( 1 : nImg ), sImg );
      x2 = x( nImg + 1 : end );
      out = wavTransH( x1 ) + curvelet( x2, 'transp' );
      out = out(:);
    end
  end

  recon = mri_reconRoemer( coilRecons );
  if strcmp( transformType, 'curvelet' ) || strcmp( transformType, 'wavCurv' )
    tmp = fdct_wrapping( recon, false );
    curvCellSizes = findCellSizes( tmp );

    if strcmp( transformType, 'wavCurv' )
      sparsifier = @(x) wavCurv( x );
      sparsifierH = @(x) wavCurv( x, 'transp' );
    else
      sparsifier = @(x) curvelet( x );
      sparsifierH = @(x) curvelet( x, 'transp' );
    end

  elseif strcmp( transformType, 'wavelet' )
    sparsifier = wavTrans;
    sparsifierH = wavTransH;

  else
    error( 'Unrecognized transform type' );
  end

  kMaskCoil = sum( abs( M ), 3 ) > 0;
  nbPerCoil = sum( kMaskCoil(:) );
  function out = A( x )
    PsiHx = reshape( sparsifierH( x ), sImg );
    FPsiHx = F( applyS( PsiHx ) );
    MFPsiHx = FPsiHx( M == 1 );
    if numel( noiseCov ) > 0
      LStarMFXin = applyL( reshape( MFPsiHx, [ nbPerCoil nCoils ] ), 'transp' );
      out = LStarMFXin(:);
    else
      out = MFPsiHx(:);
    end
  end

  MTLy = zeros( size(M) );
  function out = Aadj( y )
    if numel( noiseCov ) > 0
      Ly = applyL( reshape( y, [ nbPerCoil nCoils ] ), 'notransp' );
      MTLy( M == 1 ) = Ly;
    else
      MTLy( M == 1 ) = y;
    end
    FHMTy = applyS( FH( MTLy ), 'transp' );
    out = sparsifier( FHMTy );
  end

  applyATA = @(in) Aadj( A( in ) );

  if checkAdjoints == true
    innerProd = @(x,y) real( y(:)' * x(:) );

    xF = rand( size(kData) ) + 1i * rand( size( kData ) );
    if checkAdjoint( xF, @F, @FH, 'innerProd', innerProd ) ~= true
      error( 'FH is not the transpose of F' );
    end

    xSparsifier = rand( size(kData) ) + 1i * rand( size( kData ) );
    if checkAdjoint( xSparsifier, sparsifier, sparsifierH, 'innerProd', innerProd ) ~= true
      error( 'sparsifierH is not the transpose of sparsifier' );
    end

    xS = sparsifier( rand( sImg ) + 1i * rand( sImg ) );
    if checkAdjoint( xS, @applyS, 'innerProd', innerProd ) ~= true
      error( 'applyS adjoint test failed' );
    end

    xA = sparsifier( rand( sImg ) + 1i * rand( sImg ) );
    if checkAdjoint( xA, @A, @Aadj, 'innerProd', innerProd ) ~= true
      error( 'Aadj is not the transpose of A' );
    end
  end

  function out = g( x )
    %out = 0.5 * norm( A( x ) - beta )^2;
    out = 0.5 * norm( A( x ) - LStarBeta )^2;
  end

  AadjLStarBeta = Aadj( LStarBeta );
  function out = gGrad( x )
    out = Aadj( A( x ) ) - AadjLStarBeta;
  end

  PsiRecon0 = sparsifier( recon - reconL );  % Psi is the sparsifying transformation

  PsiRecon0 = PsiRecon0 / norm( PsiRecon0(:) ) * norm( recon(:) );
  nPsi = numel( PsiRecon0 );
  if numel( lambda ) == 0  ||  lambda == 0
    lambda = nPsi ./ ( abs( PsiRecon0 ) + reweightEpsilon );
  end

  function out = proxth( x, t )
    out = proxL1Complex( x, ( t / nPsi ) * ( w .* lambda ) );
  end

  function out = h( x )
    out = sum( abs( x(:) .* lambda(:) .* w(:) ) ) / nPsi;
  end

  if numel( stepSize ) == 0
    [normATA,piFlag] = powerIteration( applyATA, PsiRecon0, 'symmetric', true );  %#ok<ASGLU> 
    stepSize = 0.95 / normATA;
    minStep = stepSize * 1d-12;
  end

  for reweightIter = 1 : nReweightIter

    if verbose == true
      disp([ 'Working on reweighting iteration ', indx2str( reweightIter, nReweightIter ), ...
        ' of ', num2str( nReweightIter ) ]);
    end

    sMaps = mri_makeSensitivityMaps( kData, recon, 'polyOrder', polyOrder );


    if reweightIter > 1
      PsiRecon0 = xStar;
      lambda = nPsi ./ ( abs( PsiRecon0 ) + reweightEpsilon );
    end

    t = stepSize;
    if nargout > 2
      if strcmp( alg, 'pogm' )
        [xStar,oValues] = pogm( PsiRecon0, @gGrad, @proxth, nIter, 'g', @g, 'h', @h, 't', t, ...
          'tol', 1d-8, 'verbose', verbose, 'printEvery', printEvery );
      elseif strcmp( alg, 'fista' )
        [xStar,oValues] = fista( PsiRecon0, @gGrad, @proxth, 'N', nIter, 'g', @g, 'h', @h, 't', t, ...
          'tol', 1d-8, 'verbose', verbose, 'printEvery', printEvery );
      elseif strcmp( alg, 'fista_wLS' )
        [xStar,oValues] = fista_wLS( PsiRecon0, @g, @gGrad, @proxth, 'h', @h, ...
          't0', t*10, 'tol', 1d-8, 'minStep', minStep, 'N', nIter, 'verbose', verbose, ...
          'printEvery', printEvery );
      else
        error( 'Unrecognized algorithm' );
      end
    else
      if strcmp( alg, 'pogm' )
        xStar = pogm( PsiRecon0, @gGrad, @proxth, nIter, 't', t', 'tol', 1d-8 );
      elseif strcmp( alg, 'fista' )
        xStar = fista( PsiRecon0, @gGrad, @proxth, 'N', nIter, 't', t, 'tol', 1d-8 );
      elseif strcmp( alg, 'fista_wLS' )
        xStar = fista_wLS( PsiRecon0, @g, @gGrad, @proxth, 't0', t*10, 'tol', 1d-8, ...
          'minStep', minStep, 'N', nIter, 'verbose', verbose, 'printEvery', printEvery );
      else
        error( 'Unrecognized algorithm' );
      end
    end

    reconH = reshape( sparsifierH( xStar ), sImg );
    recon = reconH + reconL;
  end
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



