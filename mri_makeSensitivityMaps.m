
function senseMaps = mri_makeSensitivityMaps( kData, varargin )
  % Created using the method of "SENSE: Sensitivity Encoding for Fast MRI" by
  % Pruessmann et al., 1999
  %
  % Inputs:
  % kData - the k-space data collected from the MRI machine
  %
  % Optional Inputs:
  % L - the order of the polynomial to fit
  %
  % Outputs:
  % senseMaps - the sensitivity maps
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addParameter( 'L', 2, @ispositive );
  p.addParameter( 'mask', [], @(x) isnumeric(x) || islogical(x) );  
  p.addParameter( 'sigma', 7, @ispositive );
  p.parse( varargin{:} );
  L = p.Results.L;
  mask = p.Results.mask;
  sigma = p.Results.sigma;

  if numel( mask ) == 0, mask = ones( size( kData ) ); end
  
  coilRecons = mri_fftRecon( kData );
  ssqRecon = norms( coilRecons, 2, 4 );

  if numel( mask ) == 0, mask = ones( size( coilRecons ) ); end

  hSize = sigma * 5;
  if mod( hSize, 2 ) == 0, hSize = hSize + 1; end
  gFilt = fspecial( 'gaussian', hSize, sigma );

  coords = size2imgCoordinates( hSize );
  ys = coords(:) * ones( 1, hSize );
  xs = ones(hSize,1) * coords';

  [ nCols, nRows, nSlices, nCoils ] = size( coilRecons );

  senseMaps0 = bsxfun( @rdivide, abs( coilRecons ), ssqRecon );
  
  senseMapCols = cell( 1, nCols, 1, 1 );

  senseMaps = zeros( size( senseMaps0 ) );

  %for slice = 1 : nSlices
for slice = 10
    for i=1:nCols, senseMapCols{i} = zeros(nRows,1,nCoils); end

    parfor x0 = hSize : nCols-hSize
      disp([ 'Working on slice/col ', num2str(slice), '/', num2str(x0), ' of ', ...
        num2str(size(coilRecons,3)), '/', num2str( size( coilRecons, 2 ) ) ]);

      senseMapRowCoils = zeros( nRows, 1, nCoils );
      for y0 = hSize : nRows-hSize

        thisMask = mask( y0 - floor(hSize/2) : y0 + floor(hSize/2), ...
                         x0 - floor(hSize/2) : x0 + floor(hSize/2), slice );   %#ok<PFBNS>

        if max( thisMask(:) ) == 0, continue; end

        thisMap0 = senseMaps0( y0 - floor(hSize/2) : y0 + floor(hSize/2) , ...
                               x0 - floor(hSize/2) : x0 + floor(hSize/2), slice, : );   %#ok<PFBNS>

        thisAbsRecon = abs( ...
          coilRecons( y0 - floor(hSize/2) : y0 + floor(hSize/2) , ...
                      x0 - floor(hSize/2) : x0 + floor(hSize/2), slice, : ) ...
        );

        for coil = 1 : nCoils
          thisMap = thisMap0( :, :, :, coil );

          w = thisMask .* gFilt .* thisAbsRecon ./ abs( thisMap );
          
          c = polyFit2( xs, ys, thisMap, L, L, 'w', gFilt .* thisMask );

          senseMapRowCoils( y0, 1, coil ) = c(1,1);  % evalPoly2( c, 0, 0 );
          %senseMaps(y0,x0) = c(1,1);  % evalPoly2( c, 0, 0 );
        end
      end

      senseMapCols{x0} = senseMapRowCoils;
    end

    senseMaps(:,:,slice,:) = cell2mat( senseMapCols );
  end

  senseMaps = max( min( senseMaps, 1 ), 0 );
  
  for coil = 1 : nCoils
    senseMaps(:,:,:,coil) = mask .* senseMaps(:,:,:,coil);
  end
end
