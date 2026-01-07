
function sMaps = mrs_findSensitivityMaps( coilSpectrums, noiseCoords, varargin )
  % sMaps = mrs_findSensitivityMaps( coilSpectrums )
  %
  % Inputs:
  % coilSpectrums - an array of size sImg x nCoils x nF representing the Fourier value
  %   of each spatial location and frequency bin
  %   The data is in space and temporal frequency
  %   sImg - 2 element array representing the size of the image
  %   nCoils - the number of coils
  %   nF - the number of frequency bins in the spectrum
  % Alternatively, coilSpectrums can be sImg x nCoils x n1 x n2 x n3 x ...
  %   Each ni dimension is another dimension used for sMap estimation
  %   E.g., n1 could be nF and n2 could be number of time points
  %
  % Optional Inputs:
  % L - the order of the polynomial to use during smoothing
  % sigma - the size of the Gaussian kernel to use during smoothing
  % smoothing - (default is true) smooth the estimated sensitivity maps
  %
  % Written by Nicholas Dwork, Copyright 2021
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  p = inputParser;
  p.addParameter( 'L', 2, @ispositive );
  p.addParameter( 'sigma', 5, @ispositive );
  p.addParameter( 'smoothing', true, @(x) numel(x) == 0 || islogical(x) );
  p.parse( varargin{:} );
  L = p.Results.L;
  sigma = p.Results.sigma;
  smoothing = p.Results.smoothing;

  sImg = size( coilSpectrums, [1 2] );
  nImg = prod( sImg );
  nCoils = size( coilSpectrums, 3 );

  coilSpectrums = reshape( coilSpectrums, sImg(1), sImg(2), nCoils, [] );
  sCoilSpectrums = size( coilSpectrums );

  noise = coilSpectrums( noiseCoords(1) : noiseCoords(2), noiseCoords(3) : noiseCoords(4), :, : );  
  noise = squeeze( sqrt( sum( abs( noise ).^2, 3 ) ) );

  coilSpectrums = reshape( coilSpectrums, [ nImg sCoilSpectrums( end-1 : end ) ] );
  rsos = squeeze( sqrt( sum( abs( coilSpectrums ).^2, 2 ) ) );  % root-sum-of-squares
  denoised = removeRicianMean( rsos, noise );

  sMaps = zeros( [ nImg nCoils ] );
  for j = 1 : nImg
    %a = squeeze( sqrt( sum( coilSpectrums( j, :, : ).^2, 2 ) ) );
    a = transpose( denoised( j, : ) );
    b = transpose( squeeze( coilSpectrums( j, :, : ) ) );

    sMaps( j, : ) = a \ b;
  end

  sMaps = reshape( sMaps, [ sImg nCoils ] );

  if smoothing == true
    coilSpectrums = reshape( coilSpectrums, [ sImg nCoils, sCoilSpectrums(4) ] );
    coilRecons = sum( abs( coilSpectrums ), 4 );
    sMaps = pruessmanSmooth( sMaps, coilRecons, L, sigma );
  end

end


function smoothedMaps = pruessmanSmooth( sMaps, coilRecons, L, sigma )
  % Smoothing as in "SENSE: Sensitivity Encoding for Fast MRI" by Pruessmann et al., 1999

  [ nRows, nCols, nCoils ] = size( sMaps );

  hSize = sigma * 5;
  if mod( hSize, 2 ) == 0, hSize = hSize + 1; end
  gFilt = fspecial( 'gaussian', hSize, sigma );

  coords = size2imgCoordinates( hSize );
  ys = coords(:) * ones( 1, hSize );
  xs = ones(hSize,1) * coords';

  sMaps = padData( sMaps, [ nRows+hSize, nCols+hSize, nCoils ], 'circ', true );
  coilRecons = padData( coilRecons, [ nRows+hSize, nCols+hSize, nCoils ], 'circ', true );
  senseMapCols = cell( 1, nCols, 1 );

  %parfor x0 = floor(hSize/2)+1 : floor(hSize/2)+nCols
  for x0 = floor(hSize/2)+1 : floor(hSize/2)+nCols

    senseMapRowCoils = sMaps( :, x0, : );
    for y0 = floor(hSize/2)+1 : floor(hSize/2)+nRows

      thisMap0 = sMaps( y0 - floor(hSize/2) : y0 + floor(hSize/2) , ...
                             x0 - floor(hSize/2) : x0 + floor(hSize/2), : );   %#ok<PFBNS>

      thisRecon = coilRecons( y0 - floor(hSize/2) : y0 + floor(hSize/2) , ...
                              x0 - floor(hSize/2) : x0 + floor(hSize/2), : );   %#ok<PFBNS>

      for coil = 1 : nCoils
        thisMap = thisMap0( :, :, coil );
        thisCoilRecon = thisRecon( :, :, coil );

        %w = gFilt .* thisCoilRecon ./ abs( thisMap );
        w = gFilt .* thisCoilRecon;
        w( ~isfinite(w) ) = 0;
        c = polyFit2( xs, ys, thisMap, L, L, 'w', w );

        senseMapRowCoils( y0, 1, coil ) = c(1,1);  % evalPoly2( c, 0, 0 );
      end
    end

    senseMapCols{ x0 } = senseMapRowCoils;
  end

  smoothedMaps = cropData( cell2mat( senseMapCols ), [nRows, nCols, nCoils] );
end
