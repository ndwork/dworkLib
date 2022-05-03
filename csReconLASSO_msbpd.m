
function [ recon, oValues, lambda ] = csReconLASSO_msbpd( samples, lambda, varargin )
  % recon = csReconLASSO_msbpd( samples, lambda [, 'nIter', nIter, 'printEvery', printEvery, ...
  %   'transformType', transformType, 'verbose', verbose, 'waveletType', waveletType, ...
  %   'wavSplit', wavSplit ] )
  %
  % This routine uses proximal gradient methods to minimize 
  %   0.5 * || A x - b ||_2^2 + lambda || y ||_1
  %   where A is sampleMask * Fourier Transform * Psi.
  %   Here, Psi is either wavelets, curvelets, or both.
  %   It assumes a fully sampled center region
  %   corresponding the the two-level sampling scheme defined in
  %   "Breaking the coherence barrier: A new theory for compressed sensing"
  %   by Adcock, Ben, et al.
  %   Moreover, it uses the fully sampled scheme in a way to increase
  %   sparsity in the optimization problem according to
  %   "Utilizing the Wavelet Transform's Structure in Compressed Sensing" by Dwork et al.
  %   and "Utilizing the Structure of the Curvelet Transform with Compressed Sensing" by Dwork et al.
  %
  % Inputs:
  % samples - a 2D array that is zero wherever a sample wasn't acquired
  % lambda - regularization parameter
  %
  % Optional Inputs:
  % alg - choose the algorithm to perform the minimization.  Default is 'FISTA_wLS'.
  %   Options are FISTA_wLS.
  % nIter - the number of iterations that FISTA will perform (default is 100)
  % printEvery - FISTA prints a verbose statement every printEvery iterations
  % transformType - Choose the sparifying transformation.  Default is WavCurv (which uses
  %   the wavelet transform and the curvelet transform as a redundant dictionary).
  %   Options are:  'curvlet', 'wavelet', or 'wavCurv'.
  %     curvelet - Discrete curvelet transform
  %     wavelet - Wavelet transform (type specified by waveletType parameter)
  %     wavCurv - Redundant dictionary of wavelet and curvelet
  % verbose - if true, prints informative statements
  % waveletType - either 'Daubechies-4' for Daubechies-4 (default) or 'Haar'
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  wavSplit = zeros(4);  wavSplit(1,1) = 1;

  p = inputParser;
  p.addParameter( 'alg', 'fista_wLS', @(x) true );
  p.addParameter( 'checkAdjoints', false, @islogical );
  p.addParameter( 'nIter', [], @ispositive );
  p.addParameter( 'printEvery', 1, @ispositive );
  p.addParameter( 'stepSize', 1, @ispositive );
  p.addParameter( 'transformType', 'wavCurv', @(x) true );
  p.addParameter( 'verbose', false, @(x) isnumeric(x) || islogical(x) );
  p.addParameter( 'waveletType', [], @(x) true );
  p.addParameter( 'wavSplit', wavSplit, @isnumeric );
  p.parse( varargin{:} );
  alg = p.Results.alg;
  checkAdjoints = p.Results.checkAdjoints;
  nIter = p.Results.nIter;
  printEvery = p.Results.printEvery;
  stepSize = p.Results.stepSize;
  transformType = p.Results.transformType;
  verbose = p.Results.verbose;
  waveletType = p.Results.waveletType;
  wavSplit = p.Results.wavSplit;

  if numel( transformType ) == 0, transformType = 'wavCurv'; end
  if numel( waveletType ) == 0, waveletType = 'Daubechies-4'; end

  sImg = size( samples );
  nImg = prod( sImg );

  function out = F( x )
    out = fftshift( fftshift( ufft2( ifftshift( ifftshift( x, 1 ), 2 ) ), 1 ), 2 );
  end

  function out = FH( y )
    out = fftshift( fftshift( uifft2( ifftshift( ifftshift( y, 1 ), 2 ) ), 1 ), 2 );
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

  function [out,curvCells] = curvelet( x, type )
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

  if strcmp( transformType, 'curvelet' ) || strcmp( transformType, 'wavCurv' )
    tmp = fdct_wrapping( FH( samples ), false );
    curvCellSizes = findCellSizes( tmp );

    if strcmp( transformType, 'wavCurv' )
      sparsifier = @(x) wavCurv( x );
      sparsifierH = @(x) wavCurv( x, 'transp' );
      wavACR = makeWavAutoCalRegion( sImg, wavSplit );
      curvACR = makeCurvAutoCalRegion( sImg );
      if sum( abs( wavACR(:) ) > 0 ) > sum( abs( curvACR(:) ) > 0 )
        acr = wavACR;
      else
        acr = curvACR;
      end

    else
      sparsifier = @(x) curvelet( x );
      sparsifierH = @(x) curvelet( x, 'transp' );
      acr = makeCurvAutoCalRegion( sImg );
    end

  elseif strcmp( transformType, 'wavelet' )
    sparsifier = wavTrans;
    sparsifierH = wavTransH;
    acr = makeWavAutoCalRegion( sImg, wavSplit );

  else
    error( 'Unrecognized transform type' );
  end

  M = ( samples ~= 0 );

  function out = A( x )
    PsiHx = sparsifierH( x );
    PsiHx = reshape( PsiHx, sImg );
    FPsiHx = F( PsiHx );
    out = FPsiHx( M == 1 );
  end

  MTy = zeros( size(M) );
  function out = Aadj( y )
    MTy( M == 1 ) = y;
    FHMTy = FH( MTy );
    out = sparsifier( FHMTy );
  end

  if checkAdjoints == true
    x1 = rand( size(samples) ) + 1i * rand( size( samples ) );
    if checkAdjoint( x1, @F, @FH ) ~= true, error( 'FH is not the transpose of F' ); end

    if checkAdjoint( x1, sparsifier, sparsifierH ) ~= true
      error( 'sparsifierH is not the transpose of sparsifier' );
    end

    xA = sparsifier( double( imread( 'cameraman.png' ) ) / 255. );
    if checkAdjoint( xA, @A, @Aadj ) ~= true, error( 'Aadj is not the transpose of A' ); end
  end

  acrSamplesL = acr .* samples;
  beta = samples( M == 1 ) - acrSamplesL( M == 1 );
  function out = g( x )
    diff = A( x ) - beta;
    out = 0.5 * norm( diff(:), 2 ).^2;
  end

  AadjBeta = Aadj( beta );
  function out = gGrad( x )
    out = Aadj( A( x ) ) - AadjBeta;
  end

  %x0 = zeros( size( samples ) );
  x0 = FH( samples - acrSamplesL );
  PsiX0 = sparsifier( x0 );  % Psi is the sparsifying transformation
  nPsiX = numel( PsiX0 );
  if numel( lambda ) == 0  ||  ( lambda == 0 )
    lambda = nPsiX * findValueBelowFraction( abs( PsiX0(:) ), 0.05 );
  end

  proxth = @(x,t) proxL1Complex( x, t * lambda / nPsiX );

  function out = h( x )
    out = sum( abs( x(:) ) ) * lambda / nPsiX;
  end

  t = stepSize;
  oValues = [];
  if nargout > 1
    if strcmp( alg, 'pogm' )
      [xStar,oValues] = pogm( PsiX0, @gGrad, proxth, nIter, 'g', @g, 'h', @h, 't', t, ...
        'verbose', verbose, 'printEvery', printEvery );
    elseif strcmp( alg, 'fista' )
      [xStar,oValues] = fista( PsiX0, @gGrad, proxth, 'N', nIter, 'g', @g, 'h', @h, 't', t, ...
        'verbose', verbose, 'printEvery', printEvery );
    elseif strcmp( alg, 'fista_wLS' )
      [xStar,oValues] = fista_wLS( PsiX0, @g, @gGrad, proxth, 'h', @h, ...
        't0', t, 'N', nIter, 'verbose', verbose, 'printEvery', printEvery );
    else
      error( 'Unrecognized algorithm' );
    end
  else
    if strcmp( alg, 'pogm' )
      xStar = pogm( PsiX0, @gGrad, proxth, nIter, 't', t );
    elseif strcmp( alg, 'fista' )
      xStar = fista( PsiX0, @gGrad, proxth, 'N', nIter, 't', t );
    elseif strcmp( alg, 'fista_wLS' )
      xStar = fista_wLS( PsiX0, @g, @gGrad, proxth, 't0', t, 'N', nIter, ...
        'verbose', verbose, 'printEvery', printEvery );
    else
      error( 'Unrecognized algorithm' );
    end
  end

  blurryImg = FH( acrSamplesL );
  detailsImg = reshape( sparsifierH( xStar ), sImg );
  recon = blurryImg + detailsImg;

  showResults = false;
  if showResults == true
    cRecon = fdct_wrapping_dispcoef( fdct_wrapping( recon ) );
    cRecon( cRecon == cRecon(1,1) ) = max( cRecon(:) );
    cBlurry = fdct_wrapping_dispcoef( fdct_wrapping( blurryImg ) );
    cBlurry( cBlurry == cBlurry(1,1) ) = max( cBlurry(:) );
    cDetails = fdct_wrapping_dispcoef( fdct_wrapping( detailsImg ) );
    cDetails( cDetails == cDetails(1,1) ) = max( cDetails(:) );
    figure;  imshowscale( abs( cRecon ), 3 );  titlenice( 'cRecon' );
    figure;  imshowscale( abs( cBlurry ), 3, 'range', abs( cRecon ) );  titlenice( 'cBlurry' );
    figure;  imshowscale( abs( cDetails ), 3, 'range', abs( cRecon ) );  titlenice( 'cDetails' );

    wRecon = wavTrans( recon );
    wBlurry = wavTrans( blurryImg );
    wDetails = wavTrans( detailsImg );
    figure;  wavShow( abs( wRecon ), 3, 'wavSplit', wavSplit );  titlenice( 'wRecon' );
    figure;  wavShow( abs( wBlurry ), 3, 'wavSplit', wavSplit );  titlenice( 'wBlurry' );
    figure;  wavShow( abs( wDetails ), 3, 'wavSplit', wavSplit );  titlenice( 'wDetails' );
  end
end


function acr = makeCurvAutoCalRegion( sImg )
  [~,lowPass] = fdct_wrapping( zeros(sImg), false, 1 );
  acr = padData( lowPass, sImg );
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



