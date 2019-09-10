
function mask = mri_makeIntensityMask( kData, varargin )
  % mask = mri_makeIntensityMask( kData [, 'thresh', thresh, ...
  %   'm', m, 'morphScale', morphScale ] )
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
  p.addParameter( 'morphScale', 0, @(x) x >= 0 );
  p.addParameter( 'thresh', 0.04, @ispositive );
  p.addParameter( 'm', 2, @ispositive );
  p.parse( varargin{:} );
  morphScale = p.Results.morphScale;
  thresh = p.Results.thresh;
  m = p.Results.m;

  ssqRecon = mri_ssqRecon( kData );
  ssqRecon = ssqRecon / max( ssqRecon(:) );
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
