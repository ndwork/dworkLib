function out = fan2Parallel( sinoFan, detSize, detPlaneDist, rotDist, varargin )
  % out = fan2Parallel( sinoFan, detSize, detPlaneDist, rotDist, ...
  %                     'dProjAngle', dProjAngle, 'detCenter', detCenter, ...
  %                     'principalPoint', principalPoint, 'projAngles', projAngles )
  %
  % Inputs:
  % sinoFan - fan beam sinogram with one column per projection and one row per detector
  % detPlaneDist - the distance between the source/vertex and the detector plane
  %   Note, for optical projection tomography, this is the effective focal length
  % rotDist - the distance between the source/vertex and the center of rotation
  %
  % Optional Inputs:
  % dProjAngle - the angle between adjacent rotation angles
  %   Either dProjAngle or projAngles must be supplied
  % detCenter - detector coordinate whose ray passes through the center of rotation
  %   Defaults to the center of the detector array
  % principalPoint - detector coordinate where the principal ray intersects the detector
  %   Defaults to detCenter
  % projAngles - the angle of the principal ray for each projection
  %   Either dProjAngle or projAngles must be supplied
  %
  % Written by Nicholas Dwork, Copyright 2025
  % Modified to correctly handle principal-point offsets for real OPT data.

  p = inputParser;
  p.addParameter( 'dProjAngle', [], @ispositive );
  p.addParameter( 'detCenter', [], @isnumeric );
  p.addParameter( 'extrapolationMethod', 'none', @(x) true );
  p.addParameter( 'principalPoint', [], @isnumeric );
  p.addParameter( 'projAngles', [], @isnumeric );
  p.parse( varargin{:} );
  dProjAngle = p.Results.dProjAngle;
  detCenter = p.Results.detCenter;
  extrapolationMethod = p.Results.extrapolationMethod;
  principalPoint = p.Results.principalPoint;
  projAngles = p.Results.projAngles;

  if isempty(dProjAngle) && isempty(projAngles)
    error( 'Either dProjAngle or projAngles must be provided.' );
  end

  nDetectors = size( sinoFan, 1 );
  nProjections = size( sinoFan, 2 );

  if isempty( detCenter )
    detCenter = ( nDetectors - 1 ) / 2;
  end

  if isempty( principalPoint )
    principalPoint = detCenter;
  end

  if isempty( projAngles )
    projAngles = dProjAngle * ( 0 : nProjections - 1 );
  end
  projAngles = projAngles(:);

  % Detector coordinates measured from the principal point.
  detCoords = ( principalPoint - ( 0 : nDetectors - 1 ) ) * detSize;
  detCenterCoord = ( principalPoint - detCenter ) * detSize;

  % Fan angle of each detector ray relative to the optical axis.
  detAngles = atan2( detCoords, detPlaneDist );
  detCenterAngle = atan2( detCenterCoord, detPlaneDist );

  % Ray angle relative to the ray that passes through the center of rotation.
  rayAngles = detAngles - detCenterAngle;

  % Smallest distance from the rotation center to each ray.
  lineDists = rotDist * sin( rayAngles );

  lineDists = lineDists(:) * ones( 1, nProjections );
  lineAngles = rayAngles(:) * ones( 1, nProjections ) + ones( nDetectors, 1 ) * projAngles';

  lineAngles_ext = [ lineAngles(:) - 2*pi; lineAngles(:); lineAngles(:) + 2*pi ];
  lineDists_ext = [ lineDists(:); lineDists(:); lineDists(:) ];
  sinoFan_ext = [ sinoFan(:); sinoFan(:); sinoFan(:) ];
  F = scatteredInterpolant( lineAngles_ext, lineDists_ext, double(sinoFan_ext), 'linear', extrapolationMethod );

  maxLineDist = max( abs( lineDists(:) ) );
  nHalf = ceil( maxLineDist / detSize );
  outLineDists = ( -nHalf : nHalf ) * detSize;
  nOutLineDists = numel( outLineDists );
  outLineDists = outLineDists(:) * ones( 1, nProjections );

  outAngles = ones( nOutLineDists, 1 ) * projAngles';
  out = F( outAngles(:), outLineDists(:) );
  out = reshape( out, [ nOutLineDists, nProjections ] );
end
