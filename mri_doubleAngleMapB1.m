
function b1ScaleMap = mri_doubleAngleMapB1( dataCube, varargin )
  %b1ScaleMap = mri_mapB1( dataCube [, angles, 'mask', mask, 'verbose', verbose' ] )

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
  b1s = 0.0:0.01:1.5;
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
  simSigsSmall = cell( numel(b1s), 1 );
  simSigsLarge = cell( numel(b1s), 1 );
  %parfor b1Indx=1:numel( b1s )
for b1Indx=1:numel( b1s )
  layerSigsSmall = zeros( nVoxLayers, 1 );
  layerSigsLarge = zeros( nVoxLayers, 1 );
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
    simSigsSmall{ b1Indx } = sum( layerSigsSmall(:) );
    simSigsLarge{ b1Indx } = sum( layerSigsLarge(:) );
  end
  simSigsSmall = cell2mat( simSigsSmall );
  simSigsLarge = cell2mat( simSigsLarge );

  simDiv = abs( simSigsLarge ) ./ abs( simSigsSmall );
  b1s = b1s( isfinite( simDiv ) );
  simDiv = simDiv( isfinite( simDiv ) );

  dataDiv = abs( dataCube(:,:,2) ) ./  abs( dataCube(:,:,1) );
  b1ScaleMap = interp1( simDiv(:), b1s(:), dataDiv(:) );
  b1ScaleMap = reshape( b1ScaleMap, [ size(dataCube,1), size(dataCube,2) ] );
end

