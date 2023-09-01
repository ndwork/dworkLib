
function mask = mri_makeIntensityMask( kData, varargin )
  % mask = mri_makeIntensityMask( kData [, 'thresh', thresh, ...
  %   'm', m, 'morphScale', morphScale, 'noiseCoords', noiseCoords' ] )
  %
  % Created using the method of "SENSE: Sensitivity Encoding for Fast MRI" by
  % Pruessmann et al., 1999
  %
  % Inputs:
  % kData - the k-space data collected from the MRI machine
  %   ( kx, ky, slice, coil )
  %
  % Optional Inputs:
  % thresh - percentage of maximum intensity that is considered valid
  %   (default is 0.06)
  % m - the size of the minimum neighborhood filter
  %   (default is 5 )
  % noiseCoords - [ yMin yMax xMin xMax ] coordinates for region of image with noise
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
  p.addParameter( 'm', 2, @ispositive );
  p.addParameter( 'morphScale', 0, @(x) x >= 0 );
  p.addParameter( 'noiseCoords', [], @ispositive );
  p.addParameter( 'sigmaScalar', 3, @isnonnegative );
  p.addParameter( 'thresh', 0.04, @ispositive );
  p.parse( varargin{:} );
  m = p.Results.m;
  morphScale = p.Results.morphScale;
  noiseCoords = p.Results.noiseCoords;
  sigmaScalar = p.Results.sigmaScalar;
  thresh = p.Results.thresh;

  ssqRecon = mri_reconSSQ( kData, 'multiSlice', true );
  ssqRecon = ssqRecon / max( ssqRecon(:) );

  if numel( noiseCoords ) > 0
    noiseRegion = ssqRecon( noiseCoords(1):noiseCoords(2), noiseCoords(3):noiseCoords(4) );
    meanNoise = mean( noiseRegion(:) );
    stdNoise  = std( noiseRegion(:) );
    thresh = meanNoise + sigmaScalar * stdNoise;
  end

  mask = ssqRecon > thresh;

  if m > 0
    morphation = strel( 'disk', floor( m * morphScale ), 0 );

    nSlices = size( mask, 3 );
    maskSlices = cell( 1, 1, nSlices );
    parfor slice = 1 : nSlices

      maskSlices{slice} = ordfilt2( mask(:,:,slice), 1, ones(m) );
        % minimum neighborhood filtering

      if morphScale > 0
        maskSlices{slice} = imclose( maskSlices{slice}, morphation );

        maskSlices{slice} = imdilate( maskSlices{slice}, morphation );
      end    
    end

    mask = cell2mat( maskSlices );
  end

end
