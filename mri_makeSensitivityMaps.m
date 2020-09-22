
function senseMaps = mri_makeSensitivityMaps( kData, varargin )
  % senseMaps = mri_makeSensitivityMaps( kData [, 'L', L, ...
  %   'mask', mask, 'sigma', sigma, 'verbose', true/false ] )
  %
  % Created using the method of "SENSE: Sensitivity Encoding for Fast MRI" by
  % Pruessmann et al., 1999
  %
  % Inputs:
  % kData - the k-space data collected from the MRI machine ( kx, ky, slice, coil )
  %
  % Optional Inputs:
  % L - the order of the polynomial to fit
  % mask - an array of zeros and ones; the mask specifies those data points that have
  %  valid data (and not just noise).  It has size ( kx, ky, slice )
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

  p = inputParser;
  p.addParameter( 'L', 2, @ispositive );
  p.addParameter( 'mask', [], @(x) isnumeric(x) || islogical(x) || numel( x ) == 0 );
  p.addParameter( 'sigma', 3, @ispositive );
  p.addParameter( 'verbose', false, @(x) islogical(x) || isnumeric(x) );
  p.parse( varargin{:} );
  L = p.Results.L;
  mask = p.Results.mask;
  sigma = p.Results.sigma;
  verbose = p.Results.verbose;

  coilRecons = mri_fftRecon( kData );
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

  [ nRows, nCols, nSlices, nCoils ] = size( coilRecons );

  senseMaps0 = bsxfun( @rdivide, coilRecons, ssqRecon );
  senseMapCols = cell( 1, nCols, 1, 1 );
  senseMaps = senseMaps0;

  pfObj = parforProgress( nSlices * nCols );

  for slice = 1 : nSlices
    thisSliceMap0 = squeeze( senseMaps0(:,:,slice,:) );
    if verbose == true, disp([ 'Working on slice ', num2str(slice), ' of ', ...
      num2str(nSlices) ]); end
    for i = 1:nCols, senseMapCols{i} = thisSliceMap0(:,i,:); end

    parfor x0 = hSize : nCols-hSize
      if verbose == true, pfObj.progress( nCols*(slice-1) + x0, 20 ); end   %#ok<PFBNS>

      senseMapRowCoils = thisSliceMap0( :, x0, : );
      for y0 = hSize : nRows-hSize

        thisMask = mask( y0 - floor(hSize/2) : y0 + floor(hSize/2), ...
                         x0 - floor(hSize/2) : x0 + floor(hSize/2), slice );   %#ok<PFBNS>

        if max( thisMask(:) ) == 0, continue; end

        thisMap0 = senseMaps0( y0 - floor(hSize/2) : y0 + floor(hSize/2) , ...
                               x0 - floor(hSize/2) : x0 + floor(hSize/2), slice, : );   %#ok<PFBNS>

        thisAbsRecon = abs( ...
          coilRecons( y0 - floor(hSize/2) : y0 + floor(hSize/2) , ...
                      x0 - floor(hSize/2) : x0 + floor(hSize/2), slice, : ) );   %#ok<PFBNS>

        for coil = 1 : nCoils
          thisMap = thisMap0( :, :, :, coil );
          thisAbsCoilRecon = thisAbsRecon( :, :, :, coil );

          w = thisMask .* gFilt .* thisAbsCoilRecon ./ abs( thisMap );
          w( ~isfinite(w) ) = 0;
          c = polyFit2( xs, ys, thisMap, L, L, 'w', w );

          senseMapRowCoils( y0, 1, coil ) = c(1,1);  % evalPoly2( c, 0, 0 );
        end
      end

      senseMapCols{x0} = senseMapRowCoils;
    end

    maskedSenseMap = cell2mat( senseMapCols );
    if minMask < 1
      [ rMaskedIn,  cMaskedIn  ] = find( mask(:,:,slice) == 1 );
      [ rMaskedOut, cMaskedOut ] = find( mask(:,:,slice) == 0 );
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

    senseMaps(:,:,slice,:) = maskedSenseMap;
  end
  if verbose == true, pfObj.clean; end

  senseMaps = max( min( senseMaps, 1 ), 0 );
end
