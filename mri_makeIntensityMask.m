
function mask = mri_makeIntensityMask( kData, varargin )
  % mask = mri_makeIntensityMask( kData [, 'thresh', thresh, 'm', m ] )
  %
  % Created using the method of "SENSE: Sensitivity Encoding for Fast MRI" by
  % Pruessmann et al., 1999
  %
  % Inputs:
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
  p.addParameter( 'thresh', 0.06, @ispositive );
  p.addParameter( 'm', 5, @ispositive );
  p.parse( varargin{:} );
  thresh = p.Results.thresh;
  m = p.Results.m;

  if nargin < 2, thresh = 0.06; end

  ssqRecon = mri_ssqRecon( kData );
  ssqRecon = ssqRecon / max( ssqRecon(:) );
  mask = ssqRecon > thresh;

  for slice = 1 : size( mask, 3 )
    mask(:,:,slice) = ordfilt2( mask(:,:,slice), 1, ones(m) );
      % minimum neighborhood filtering
  end
end
