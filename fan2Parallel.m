
function out = fan2Parallel( sinoFan, detSize, detPlaneDist, rotDist, varargin )
  % out = fan2Parallel( sinoFan, focalLength, rotDist )
  %
  % Inputs:
  % sinoFan - fan beam sinogram with one column per projection and one row per detector
  % detPlaneDistb - the distance between the vertex and the detector
  %   Note, for optical projection tomography, this is the focal length
  % rotDist - the distance between the vertex and the center of rotation
  %
  % Optional Inputs:
  % dProjAngle - the angle between adjacent rotation angles
  %   Either dAngle or projAngles must be supplied
  % detCenter - the location of the center detector (defaults to center of projection)
  %   Note, for optical projection tomography, this is a component of the principal point
  % projAngles - the angle of the principal ray for each projection
  %   Either dAngle or projAngles must be supplied
  %
  % Written by Nicholas Dwork, Copyright 2025
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addParameter( 'dProjAngle', [], @ispositive );
  p.addParameter( 'detCenter', [], @isnumeric );
  p.addParameter( 'extrapolationMethod', 'none', @(x) true );
  p.addParameter( 'projAngles', [], @isnumeric );
  p.parse( varargin{:} );
  dProjAngle = p.Results.dProjAngle;
  detCenter = p.Results.detCenter;
  extrapolationMethod = p.Results.extrapolationMethod;
  projAngles = p.Results.projAngles;

  if isempty(dProjAngle) && isempty(projAngles)
    error( 'Either dProjAngle or projAngles must be provided.' );
  end

  nDetectors = size( sinoFan, 1 );
  nProjections = size( sinoFan, 2 );

  if numel( detCenter ) == 0
    detCenter = ( nDetectors - 1 ) / 2;
  end

  if numel( projAngles ) == 0
    projAngles = dProjAngle * ( 0 : nProjections -1 );
  end

  detLocs = -( ( 0 : nDetectors-1 ) - detCenter ) * detSize;
  detAngles = atan2( detLocs, detPlaneDist );
  lineDists = rotDist * sin( detAngles );   % Smallest distance from rotation center to projection line

  lineDists = lineDists(:) * ones(1,nProjections );
  lineAngles = detAngles(:) * ones(1,nProjections) + ones(nDetectors,1) * projAngles(:)';

  lineAngles_ext = [lineAngles(:) - 2*pi; lineAngles(:); lineAngles(:) + 2*pi ];
  lineDists_ext = [lineDists(:); lineDists(:); lineDists(:)];
  sinoFan_ext = [sinoFan(:); sinoFan(:); sinoFan(:)];
  F = scatteredInterpolant(lineAngles_ext, lineDists_ext, sinoFan_ext, 'linear', extrapolationMethod );

  maxLineDist = max( abs( lineDists(:) ) );
  nOutLineDists = ceil( 2 * maxLineDist / detSize );
  outLineDists = linspace( -maxLineDist, maxLineDist, nOutLineDists );
  outLineDists = outLineDists(:) * ones(1, nProjections);

  outAngles = ones(nOutLineDists,1) * projAngles(:)';
  out = F( outAngles(:), outLineDists(:) );
  out = reshape( out, [ nOutLineDists nProjections ] );
end
