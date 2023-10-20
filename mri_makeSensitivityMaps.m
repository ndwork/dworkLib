
function senseMaps = mri_makeSensitivityMaps( kData, varargin )
  % Either
  %   senseMaps = mri_makeSensitivityMaps( kData, 'alg', 'pruessman', [, 'L', L, ...
  %     'mask', mask, 'sigma', sigma, 'verbose', true/false ] )
  % or
  %   senseMaps = mri_makeSensitivityMaps( kData, 'alg', 'simple' [, 'epsilon', epsilon ] );
  % or
  %   senseMaps = mri_makeSensitivityMaps( kData, img, 'alg', 'ying' [, 'polyOrder', polyOrder, ...
  %     'verbose', true/false );
  %
  % Options are to either use the method of (1) or (2, default)
  % 1) Created using the method of "SENSE: Sensitivity Encoding for Fast MRI" by
  %    Pruessmann et al., 1999
  % 2) Created using the method of "Joint image reconstruction and sensitivity estimation "
  %    in SENSE (JSENSE)" by Ying et al.
  %
  % Inputs:
  % kData - the k-space data collected from the MRI machine ( kx, ky, coil )
  %
  % Optional Inputs:
  % epsilon - a term added to the denominator to prevent divide by 0
  % L - the order of the polynomial to fit
  % mask - an array of zeros and ones; the mask specifies those data points that have
  %  valid data (and not just noise).  It has size ( kx, ky )
  % sigma - the standard deviation of the Gaussian weights to use when fitting a polynomial
  % verbose - if set to true, displays informative statements
  %
  % Outputs:
  % senseMaps - the sensitivity maps
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp([ 'Usage:  senseMaps = mri_makeSensitivityMaps( kData [, ''L'', L, ''mask'', mask, ', ...
      '''sigma'', sigma, ''verbose'', true/false ] ) ' ]);
    if nargout > 0, senseMaps=[]; end
    return
  end

  defaultPolyOrder = 17;

  p = inputParser;
  p.addOptional( 'img', [] );
  p.addParameter( 'alg', 'ying', @(x) true );
  p.addParameter( 'epsilon', 0, @ispositive );
  p.addParameter( 'L', 2, @ispositive );
  p.addParameter( 'mask', [], @(x) isnumeric(x) || islogical(x) || numel( x ) == 0 );
  p.addParameter( 'polyOrder', defaultPolyOrder, @ispositive );
  p.addParameter( 'sigma', 3, @ispositive );
  p.addParameter( 'verbose', false, @(x) islogical(x) || isnumeric(x) );
  p.parse( varargin{:} );
  alg = p.Results.alg;
  epsilon = p.Results.epsilon;
  img = p.Results.img;
  L = p.Results.L;
  mask = p.Results.mask;
  polyOrder = p.Results.polyOrder;
  sigma = p.Results.sigma;
  verbose = p.Results.verbose;

  if numel( epsilon ) == 0, epsilon = 0; end
  if numel( polyOrder ) == 0, polyOrder = defaultPolyOrder; end

  switch alg
    case 'pruessman'
      senseMaps = mri_makeSensitivityMaps_pruessman( kData, L, mask, sigma, verbose );
    case 'simple'
      senseMaps = mri_makeSensitivityMaps_simple( kData, epsilon );
    case 'sortaSimple'
      senseMaps = mri_makeSensitivityMaps_sortaSimple( kData, epsilon );
    case 'ying'
      senseMaps = mri_makeSensitivityMaps_ying( kData, img, polyOrder );
    otherwise
      error( 'Unrecognized algorithm' );
  end

end

function senseMaps = mri_makeSensitivityMaps_pruessman( kData, L, mask, sigma, verbose )

  coilRecons = mri_reconIFFT( kData, 'multiSlice', true );
  sCoilRecons = size( coilRecons );
  ssqRecon = norms( coilRecons, 2, 4 );

  if numel( mask ) == 0
    mask = ones( sCoilRecons(1:3) );
    minMask = 1;
  else
    minMask = min( mask(:) );
  end

  hSize = sigma * 5;
  if mod( hSize, 2 ) == 0, hSize = hSize + 1; end
  gFilt = fspecial( 'gaussian', hSize, sigma );

  coords = size2imgCoordinates( hSize );
  ys = coords(:) * ones( 1, hSize );
  xs = ones(hSize,1) * coords';

  [ nRows, nCols, nCoils ] = size( coilRecons );

  senseMaps0 = bsxfun( @rdivide, coilRecons, ssqRecon );
  senseMaps = zeros( size( senseMaps0 ) );

  senseMapCols = cell( 1, nCols, 1 );
  for i = 1:nCols, senseMapCols{i} = senseMaps0(:,i,:); end

  pfObj = parforProgress( nCols );
  parfor x0 = hSize : nCols-hSize
    if verbose == true, pfObj.progress( nCols + x0, 20 ); end   %#ok<PFBNS>

    senseMapRowCoils = senseMaps0( :, x0, : );
    for y0 = hSize : nRows-hSize

      thisMask = mask( y0 - floor(hSize/2) : y0 + floor(hSize/2), ...
                       x0 - floor(hSize/2) : x0 + floor(hSize/2) );   %#ok<PFBNS>

      if max( thisMask(:) ) == 0, continue; end

      thisMap0 = senseMaps0( y0 - floor(hSize/2) : y0 + floor(hSize/2) , ...
                               x0 - floor(hSize/2) : x0 + floor(hSize/2), : );   %#ok<PFBNS>

      thisAbsRecon = abs( ...
        coilRecons( y0 - floor(hSize/2) : y0 + floor(hSize/2) , ...
                    x0 - floor(hSize/2) : x0 + floor(hSize/2), : ) );   %#ok<PFBNS>

      for coil = 1 : nCoils
        thisMap = thisMap0( :, :, coil );
        thisAbsCoilRecon = thisAbsRecon( :, :, coil );

        w = thisMask .* gFilt .* thisAbsCoilRecon ./ abs( thisMap );
        w( ~isfinite(w) ) = 0;
        c = polyFit2( xs, ys, thisMap, L, L, 'w', w );

        senseMapRowCoils( y0, 1, coil ) = c(1,1);  % evalPoly2( c, 0, 0 );
      end
    end

    senseMapCols{ x0 } = senseMapRowCoils;
  end

  maskedSenseMap = cell2mat( senseMapCols );
  if minMask < 1
    [ rMaskedIn,  cMaskedIn  ] = find( mask(:,:) == 1 );
    [ rMaskedOut, cMaskedOut ] = find( mask(:,:) == 0 );
    nMaskedIn = numel( rMaskedIn );
    nMaskedOut = numel( rMaskedOut );

    for coilIndx = 1 : nCoils
      sMaskedSenseMap = size(maskedSenseMap);
      theseIndxs = sub2ind( sMaskedSenseMap, rMaskedIn, cMaskedIn, ...
        coilIndx*ones(nMaskedIn,1) );
      theseSensed = maskedSenseMap( theseIndxs );
      SI = scatteredInterpolant( cMaskedIn, rMaskedIn, theseSensed, 'linear', 'nearest' );

      thoseIndxs = sub2ind( sMaskedSenseMap, rMaskedOut, cMaskedOut, coilIndx*ones(nMaskedOut,1) );
      maskedSenseMap( thoseIndxs ) = SI( cMaskedOut, rMaskedOut );
    end
  end

  senseMaps(:,:,:) = maskedSenseMap;

  if verbose == true, pfObj.clean; end
end


function senseMaps = mri_makeSensitivityMaps_simple( kData, epsilon )
  coilRecons = mri_reconIFFT( kData, 'multiSlice', true );
  reconSSQ = sqrt( sum( coilRecons .* conj( coilRecons ), ndims( kData ) ) );
  senseMaps = bsxfun( @rdivide, coilRecons, ( reconSSQ + epsilon ) );
end


function sMaps = mri_makeSensitivityMaps_sortaSimple( kData, epsilon )

  nCoils = size( kData, 3 );
  mask = mri_makeIntensityMask( kData, 'noiseCoords', [1 25 1 25], 'sigmaScalar', 5 );
  if max( mask(:) ) == 0
    sMaps = zeros( kData );
    return
  end
  sMaps = mri_makeSensitivityMaps_simple( kData, epsilon );
  sMaps = bsxfun( @times, sMaps, mask );
  
  [ rMaskedIn,  cMaskedIn  ] = find( mask == 1 );
  [ rMaskedOut, cMaskedOut ] = find( mask == 0 );
  sMask = size( mask );
  theseIndxs = sub2ind( sMask, rMaskedIn, cMaskedIn );
  thoseIndxs = sub2ind( sMask, rMaskedOut, cMaskedOut );

  for coilIndx = 1 : nCoils
    maskedSenseMap = sMaps( :, :, coilIndx );
    theseSensed = maskedSenseMap( theseIndxs );
    SI = scatteredInterpolant( cMaskedIn, rMaskedIn, theseSensed, 'natural', 'nearest' );     
    maskedSenseMap( thoseIndxs ) = SI( cMaskedOut, rMaskedOut );
    sMap = smoothImg( maskedSenseMap, 'gaussian', 5 );
    sMaps( :, :, coilIndx ) = sMap;
  end
end


function sMaps = mri_makeSensitivityMaps_ying( kData, img, polyOrder )

  if numel( polyOrder ) == 1, polyOrder = [ polyOrder polyOrder ]; end

  sKData = size( kData );
  nPixels = prod( sKData(1:2) );
  nCoils = sKData( 3 );
  coords = size2imgCoordinates( sKData(1:2) );
  [ ys, xs ] = ndgrid( coords{1} / max(coords{1}), coords{2} / max( coords{2} ) );

  kMask = max( abs( kData ), [], 3 ) ~= 0;

  % Make the 2D Vandermonde matrix
  V = zeros( sKData(1), sKData(2), prod( polyOrder+1 ) );
  orderIndx = 0;
  for xOrder = 0 : polyOrder(2)
    for yOrder = 0 : polyOrder(1)
      orderIndx = orderIndx + 1;
      V(:,:,orderIndx) = ( xs.^xOrder ) .* ( ys.^yOrder );
    end
  end
  Vmat = reshape( V, nPixels, [] );

  function out = applyA( in, type )
    if nargin < 2 || strcmp( type, 'notransp' )
      in = reshape( in, size( V, 3 ), [] );
      Vin = reshape( Vmat * in, sKData );
      imgVin = bsxfun( @times, Vin, img );
      out = fftshift2( fft2( ifftshift2( imgVin ) ) );
      out = reshape( out, nPixels, [] );
      out = out( kMask ~= 0, : );
    else
      kIn = zeros( sKData );
      kIn( abs( kData ) ~= 0 ) = in;
      FhkIn = fftshift2( fft2h( ifftshift2( kIn ) ) );
      imgFhkIn = bsxfun( @times, FhkIn, conj(img) );
      out = Vmat' * reshape( imgFhkIn, nPixels, [] );
    end
    out = out(:);
  end

  doCheckAdjoint = false;
  if doCheckAdjoint == true
    [checkA,errA] = checkAdjoint( rand(324,8), @applyA );
    if checkA ~= true, error([ 'Adjoint check failed with error: ', num2str(errA) ]); end
  end

  polyCoeffs = lsqr( @applyA, kData( abs(kData) ~= 0 ), [], 100, [] );
  polyCoeffs = reshape( polyCoeffs, [], nCoils );

  sMaps = tensorprod( V, polyCoeffs, 3, 1 );
end



