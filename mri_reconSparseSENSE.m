
function recon = mri_reconSparseSENSE( kData, sMaps, lambda, varargin )
  % recon = mri_reconSparseSENSE( kData, sMaps[, 'img0', img0, 'nIter', nIter, ...
  %   'noiseCov', noiseCov ] )
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
  p.addRequired( 'kData', @isnumeric );
  p.addRequired( 'sMaps', @isnumeric );
  p.addRequired( 'lambda', @ispositive );
  p.addParameter( 'checkAdjoints', false, @islogical );
  p.addParameter( 'debug', false, @(x) isnumeric(x) || islogical(x) );
  p.addParameter( 'img0', [], @isnumeric );
  p.addParameter( 'nIter', [], @ispositive );
  p.addParameter( 'noiseCov', [], @isnumeric );
  p.addParameter( 'optAlg', 'fista_wLS', @(x) true );
  p.addParameter( 'printEvery', 10, @ispositive );
  p.addParameter( 'verbose', false, @(x) isnumeric(x) || islogical(x) );
  p.addParameter( 'waveletType', 'Daubechies', @(x) true );
  p.parse( kData, sMaps, lambda, varargin{:} );
  checkAdjoints = p.Results.checkAdjoints;
  debug = p.Results.debug;
  img0 = p.Results.img0;
  nIter = p.Results.nIter;
  noiseCov = p.Results.noiseCov;
  optAlg = p.Results.optAlg;
  printEvery = p.Results.printEvery;
  waveletType = p.Results.waveletType;
  verbose = p.Results.verbose;

  if numel( nIter ) == 0
    if debug == true
      nIter = 30;
    else
      nIter = 100;
    end
  end

  function out = applyS( in, type )
    if nargin < 2 || strcmp( type, 'notransp' )
      out = bsxfun( @times, sMaps, in );
    else
      out = sum( conj(sMaps) .* in, 3 );
    end
  end

  function out = applyF( in, type )
    if nargin < 2 || strcmp( type, 'notransp' )
      out = fftshift2( fft2( ifftshift2( in ) ) );
    else
      out = fftshift2( fft2h( ifftshift2( in ) ) );
    end
  end

  function out = applySF( in, type )
    if nargin < 2 || strcmp( type, 'notransp' )
      out = applyF( applyS( in ) );
    else
      out = applyS( applyF( in, 'transp' ), 'transp' );
    end
  end

  sKData = size( kData );
  dataMask = ( abs(kData) ~= 0 );
  function out = applyA( in, type )
    if nargin < 2 || strcmp( type, 'notransp' )
      in = reshape( in, sKData(1:2) );
      out = applyF( applyS( in ) );
      out = out( dataMask == 1 );
    else
      tmp = zeros( sKData );
      tmp( dataMask == 1 ) = in;
      in = tmp;  clear tmp;
      in = reshape( in, sKData );
      out = applyS( applyF( in .* dataMask, 'transp' ), 'transp' );
    end
  end

  applyATA = @(x) applyA( applyA( x ), 'transp' );

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

  M = max( dataMask, [], 3 );
  nM = sum( M(:) );
  nCoils = sKData( ndims( kData ) );
  b = zeros( nM, nCoils );
  for coilIndx = 1 : nCoils
    coilSamples = kData( :, :, coilIndx );
    b( :, coilIndx ) = coilSamples( M == 1 );
  end
  b = b(:);
 
  function out = g( x )
    diff = applyA( x ) - b;
    if numel( noiseCov ) > 0
      NInvDiff = applyL( reshape( diff, [ nM nCoils ] ) );
      out = real( 0.5 * diff' * NInvDiff(:) );
    else
      out = 0.5 * norm( diff(:), 2 )^2;
    end
  end

  if numel( noiseCov ) > 0
    NInvb = applyL( reshape( b, [ nM nCoils ] ) );
    AadjMb = applyA( NInvb(:), 'transp' );
  else
    Aadjb = applyA( b, 'transp' );
  end
  function out = gGrad( x )
    if numel( noiseCov ) > 0
      Ax = applyA( x );
      NInvAx = applyL( reshape( Ax, [ nM nCoils ] ) );
      out = applyA( NInvAx(:), 'transp' ) - AadjMb;
    else
      out = applyA( applyA( x ), 'transp' ) - Aadjb;
    end
  end


  if numel( img0 ) == 0
    coilRecons = mri_reconIFFT( kData );
    img0 = mri_reconRoemer( coilRecons );
  end
  split = makeWavSplit( size( img0 ) );

  if strcmp( waveletType, 'Daubechies' )
    wavOp = @(x) wtDaubechies2( x, split );
    wavAdj = @(y) iwtDaubechies2( y, split );

  elseif strcmp( waveletType, 'Haar' )
    wavOp = @(x) wtHaar2( x, split );
    wavAdj = @(y) iwtHaar2( y, split );

  else
    error( 'Unrecognized wavelet type' );
  end

  if checkAdjoints == true
    innerProd = @(x,y) real( dotP( x, y ) );
    [checkS,errS] = checkAdjoint( img0, @applyS, 'innerProd', innerProd );
    if checkS ~= 1, error( ['Adjoint of S failed with error ', num2str(errS) ]); end
    [checkF,errF] = checkAdjoint( repmat(img0,[1,1,8]), @applyF, 'innerProd', innerProd );
    if checkF ~= 1, error( ['Adjoint of F failed with error ', num2str(errF) ]); end
    [checkSF,errSF] = checkAdjoint( img0, @applySF, 'innerProd', innerProd );
    if checkSF ~= 1, error( ['Adjoint of SF failed with error ', num2str(errSF) ]); end
    [checkA,errA] = checkAdjoint( img0, @applyA, 'innerProd', innerProd );
    if checkA ~= 1, error( ['Adjoint of A failed with error ', num2str(errA) ]); end
    if checkAdjoint( imgRand, wavOp, wavAdj ) ~= true, error( 'WT is not the transpose of W' ); end
    if checkW ~= 1, error( ['Adjoint of wavOp failed with error ', num2str(errW) ]); end
  end

  function out = proxth( x, t )
    out = wavAdj( softThresh( wavOp(x), t*lambda ) );
  end

  function out = h( x )
    Wx = wavOp(x);
    out = lambda * sum( abs( Wx(:) ) );
  end

  normATA = powerIteration( applyATA, rand( size( img0 ) ), 'symmetric', true );

  if normATA == 0
    % A just sets everything to 0
    recon = zeros( size( img0 ) );
    oValues = g( recon ) * ones( nIter, 1 );   %#ok<NASGU>
    return;
  end

  t = 0.99 / normATA;

  if debug
    if strcmp( optAlg, 'fista' )
      [recon,oValues,relDiffs] = fista( img0, @gGrad, @proxth, 't', t, ...
        'g', @g, 'h', @h, 'printEvery', printEvery, 'verbose', verbose );   %#ok<ASGLU>
    elseif strcmp( optAlg, 'fista_wLS' )
      t = t * 10;
      [recon,oValues] = fista_wLS( img0, @g, @gGrad, @proxth, 'h', @h, ...
        't0', t, 'N', nIter, 'restart', true, 'verbose', true, 'printEvery', printEvery );   %#ok<ASGLU>
    elseif strcmp( optAlg, 'pogm' )
      [recon,oValues] = pogm( img0, @gGrad, @proxth, nIter, 't', t, 'g', @g, 'h', @h, ...
        'printEvery', printEvery, 'verbose', true );   %#ok<ASGLU>
    else
      error([ 'Unrecognized optAlg: ', optAlg ]);
    end
  else
    if strcmp( optAlg, 'fista' )
      recon = fista( img0, @gGrad, @proxth, 't', t, 'printEvery', printEvery, 'verbose', true );
    elseif strcmp( optAlg, 'fista_wLS' )
      t = t * 10;
      recon = fista_wLS( img0, @g, @gGrad, @proxth, 't0', t, 'N', nIter, ...
        'restart', true, 'verbose', verbose, 'printEvery', printEvery );
    elseif strcmp( optAlg, 'pogm' )
      recon = pogm( img0, @gGrad, @proxth, nIter, 't', t, 'g', @g, 'h', @h, ...
        'printEvery', printEvery, 'verbose', verbose );
    else
      error([ 'Unrecognized optAlg: ', optAlg ]);
    end
  end
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


