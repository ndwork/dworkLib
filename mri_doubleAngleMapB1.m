
function b1ScaleMap = mri_doubleAngleMapB1( dataCube, sliceThickness, varargin )
  % b1ScaleMap = mri_doubleAngleMapB1( dataCube, sliceThickness [, angles, 'mask', mask, 'alg', alg, ...
  %   'verbose', verbose' ] )
  %
  % Performs double angle B1 mapping
  %
  % Inputs:
  % dataCube - a 3D array of size MxNx2.  The first / second slice is the
  %   single / double angle data.
  % sliceThickness - used when designing RF pulse (in cm)
  %
  % Optional Inputs:
  % angles - two element array specifying the flip angles of the data (radians)
  % mask - a 2D array of size MxN.  Only map pixels with nonzero mask values.
  % simple - Either 0 or 1
  %   If 0, fits an exponential decay to the data
  %   If 1, takes the slice profile into account
  % verbose - scalar; info statements made if verbose is nonzero
  %
  % Outputs:
  % b1ScaleMap - a 2D array of size MxN; values are b1 scaling
  %
  % Written by Nicholas Dwork - Copyright 2018
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.


  p = inputParser;
  p.addOptional( 'angles', [pi/3 2*pi/3], @isnumeric );
  p.addParameter( 'mask', [], @(x) isnumeric(x) || islogical(x) );
  p.addParameter( 'alg', [], @(x) true );
  p.addParameter( 'verbose', 0, @(x) isnumeric(x) || islogical(x) );
  p.parse( varargin{:} );
  angles = p.Results.angles;
  mask = p.Results.mask;
  alg = p.Results.alg;
  verbose = p.Results.verbose;

  if size( dataCube, 3 ) ~= 2, error('mri_doubleAngleMapB1: must supply two images'); end
  if numel( angles ) ~= 2, error('mri_doubleAngleMapB1: must supply two angles'); end

  if strcmp( alg, 'simple' )

    singleAngleImg = abs( dataCube(:,:,1) );
    doubleAngleImg = abs( dataCube(:,:,2) );
    tmp = doubleAngleImg ./ ( 2 * singleAngleImg );
    angleMap = acos( max( min( tmp, 1 ), -1 ) );
    b1ScaleMap = angleMap / angles(1);

  else

    tbw = 8;
    rfDuration = 1.28;  % ms
    b1s = 0.01:0.01:1.5;
    B0 = 0;  % residual B0 offset

    nRF = tbw * 100 + 1;
    rf = real( dzrf( nRF, tbw, 'ex', 'ms' ) )';  % Hamming windowed sinc function
    rfSmall = rf / sum( rf ) * angles(1);
    rfLarge = rf / sum( rf ) * angles(2);

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
    parfor b1Indx=1:numel( b1s )
      layerSigsSmall = zeros( nVoxLayers, 1 );
      layerSigsLarge = zeros( nVoxLayers, 1 );
      if verbose ~= 0 && mod( b1Indx, 10 ) == 0
        disp([ 'Working on ', num2str(b1Indx), ' of ', numel( b1s ) ]);
      end
      b1 = b1s( b1Indx );
      %[aSmall,bSmall] = abr( b1*rfSmall_total, gz_total, zSlice, B0 );
      %[aLarge,bLarge] = abr( b1*rfLarge_total, gz_total, zSlice, B0 );
      [aSmall,bSmall] = abrm( b1*rfSmall_total, gz_total, zSlice, B0 );
      [aLarge,bLarge] = abrm( b1*rfLarge_total, gz_total, zSlice, B0 );
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
    b1ScaleMap( ~isfinite( b1ScaleMap ) ) = 0;
    b1ScaleMap = reshape( b1ScaleMap, [ size(dataCube,1), size(dataCube,2) ] );
    
  end

  if numel( mask ) > 0,  b1ScaleMap = b1ScaleMap .* mask; end
end

