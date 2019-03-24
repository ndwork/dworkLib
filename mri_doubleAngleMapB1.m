
function angleMap = mri_doubleAngleMapB1( dataCube, varargin )
  %b1Map = mri_mapB1( dataCube [, angles, 'mask', mask, 'verbose', verbose' ] )

  p = inputParser;
  p.addOptional( 'angles', [60 120], @isnumeric );
  p.addParameter( 'mask', [], @(x) isnumeric(x) || islogical(x) );
  p.addParameter( 'verbose', 0, @(x) isnumeric(x) || islogical(x) );
  p.parse( varargin{:} );
  angles = p.Results.angles;
  mask = p.Results.mask;
  verbose = p.Results.verbose;

  if size( dataCube, 3 ) ~= 2, error('mri_mapB1: must supply two images'); end
  if numel( angles ) ~= 2, error('mri_mapB1: must supply two angles'); end
  if numel( mask ) == 0, mask = ones( size(dataCube(:,:,1)) ); end

  tbw = 8;
  rfDuration = 1.28;  % ms
  sliceThickness = 0.5; % cm
  b1s = 0.0:0.01:180/angles(2);
  B0 = 0;  % residual B0 offset

  nRF = tbw * 100 + 1;
  rf = real( dzrf( nRF, tbw, 'ex', 'ms' ) )';  % Hamming windowed sinc function
  rfLarge = rf / sum( rf ) * angles(2) * pi/180;
  rfSmall = rf / sum( rf ) * angles(1) * pi/180;

  dtRF = rfDuration / (nRF-1);  % dtRf is time between samples of pulse
  nVoxLayers = tbw * 10 + 1;  % Number of layers to break the voxel up into
  zL = -sliceThickness * (4/2);  % Rule of thumb: simulate for 4 x slice thickness
  zH =  sliceThickness * (4/2);
  layerIndxs = -floor(nVoxLayers/2) : floor((nVoxLayers-1)/2);
  dz = ( zH - zL ) / ( nVoxLayers + mod(nVoxLayers+1,2) );
  zSlice = dz * layerIndxs;

  [gammaBar,gamma] = getGammaH();  % kHz/Gauss, kRad/Gauss
  bw = tbw / rfDuration;  % kHz
  gzrf_amplitude = bw / gammaBar / sliceThickness;  % Gauss per cm
  gz_rot = gamma * gzrf_amplitude * dtRF * ones( nRF, 1 );	% radians/cm
  splitIndx = round( numel(gz_rot)/2 );
  gzRefocus = -gz_rot( splitIndx:end );
  gz_total = [ gz_rot; gzRefocus; ];
  rfSmall_total = [ rfSmall(:); zeros(size(gzRefocus)); ];
  rfLarge_total = [ rfLarge(:); zeros(size(gzRefocus)); ];

  R2C = mri_makeSigConversionMatrix_r2c();
  C2R = mri_makeSigConversionMatrix_c2r();
  layerSigsSmall = zeros( nVoxLayers, 1 );
  layerSigsLarge = zeros( nVoxLayers, 1 );
  sigsSmall = zeros( numel(b1s), 1 );
  sigsLarge = zeros( numel(b1s), 1 );
  for b1Indx=1:numel( b1s )
    if verbose ~= 0 && mod( b1Indx, 10 ) == 0
      disp([ 'Working on ', num2str(b1Indx), ' of ', numel( b1s ) ]);
    end
    b1 = b1s( b1Indx );
    [aSmall,bSmall] = abr( b1*rfSmall_total, gz_total, zSlice, B0 );
    [aLarge,bLarge] = abr( b1*rfLarge_total, gz_total, zSlice, B0 );
    rfMsSmall = mri_ab2Matrix( aSmall, bSmall );
    rfMsLarge = mri_ab2Matrix( aLarge, bLarge );
    for layer=1:nVoxLayers
      MprimeSmall = C2R * rfMsSmall(:,:,layer) * R2C * [0; 0; 1];
      MprimeLarge = C2R * rfMsLarge(:,:,layer) * R2C * [0; 0; 1];
      layerSigsSmall( layer ) = MprimeSmall(1) + 1i * MprimeSmall(2);
      layerSigsLarge( layer ) = MprimeLarge(1) + 1i * MprimeLarge(2);
    end
    sigsSmall( b1Indx ) = mean( layerSigsSmall );
    sigsLarge( b1Indx ) = mean( layerSigsLarge );
  end
  divDataSim = abs( sigsLarge ) ./ abs( sigsSmall );
  simAngles = acos( 0.5 * divDataSim );
  simAngles( ~isfinite( simAngles ) ) = 0;

  idealSigsSmall = sin( angles(1) * pi/180 * b1s );
  idealSigsLarge = sin( angles(2) * pi/180 * b1s );
  idealAngles = acos( 0.5 * idealSigsLarge ./ idealSigsSmall );
  idealAngles( ~isfinite(idealAngles) ) = 0;

  divData = abs(dataCube(:,:,2)) ./ abs(dataCube(:,:,1)) .* mask;
  divData = max( min( divData, 2 ), 0 );
  angleMapIdeal = acos( 0.5 * divData ) .* mask;

  angleMap = interp1( idealAngles, simAngles, angleMapIdeal(:) );
  angleMap = reshape( angleMap, size( angleMapIdeal ) );
  %figure; imshowscale( angleMap*180/pi, 3 ); titlenice( 'Angle Map with slice profile' );

disp('NICK, YOU STILL NEED TO GET THE B1 MAP FROM THE ANGLE MAP');
  
  %error('How do I get the B1 map from the angle map?');
end

