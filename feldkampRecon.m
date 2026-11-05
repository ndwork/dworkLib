function out = feldkampRecon( sinoFlatPanel, sImg, f, d2Rot, varargin )
  % Approximate reconstruction for 3D computed tomography (CT) using the feldkamp method
  %
  % Inputs:
  % sinoFlatPanel - a 3D array of size nY x nX x nProjections representing data collected by a flat panel
  % sImg - a 2D array specifying the size of the output image for each slice of the reconstruction
  % f - the distance from the fulcrum to the detector plane
  % d2Rot - the distance from the fulcrum to the center of rotation
  %
  % Optional Inputs:
  % dProjAngle - the angle between adjacent rotation angles
  %   Either dProjAngle or projAngles must be supplied
  % dSizeX - the size of the detector in the dimension of rotation
  % dSizeY - the size of the detector in the dimension of the rotation axis
  % detCenter - detector coordinate whose ray passes through the center of rotation
  % inputIsLineIntegrals - if true, sinoFlatPanel is already -log(I/I0)
  % I0 - scalar or array used when inputIsLineIntegrals is false
  % principalPoint - detector coordinate where the principal ray intersects the detector
  % projAngles - the angle of the principal ray for each projection
  %   Either dProjAngle or projAngles must be supplied
  % verbose - true/false to print verbosity statements
  %
  % Written by Nicholas Dwork, Copyright 2025
  % Modified to propagate real-data geometry and accept pre-logged line integrals.

  p = inputParser;
  p.addParameter( 'cropAmount', [], @ispositive );
  p.addParameter( 'dProjAngle', [], @isnumeric );
  p.addParameter( 'dSizeX', 1, @ispositive );
  p.addParameter( 'dSizeY', 1, @ispositive );
  p.addParameter( 'detCenter', [], @isnumeric );
  p.addParameter( 'I0', 1, @isnumeric );
  p.addParameter( 'inputIsLineIntegrals', false, @islogical );
  p.addParameter( 'outDir', [], @(x) true );
  p.addParameter( 'printEvery', 10, @ispositive );
  p.addParameter( 'principalPoint', [], @isnumeric );
  p.addParameter( 'projAngles', [], @isnumeric );
  p.addParameter( 'verbose', false, @islogical );
  p.parse( varargin{:} );
  cropAmount = p.Results.cropAmount;
  dProjAngle = p.Results.dProjAngle;
  dSizeX = p.Results.dSizeX;
  dSizeY = p.Results.dSizeY; %#ok<NASGU>
  detCenter = p.Results.detCenter;
  I0 = p.Results.I0;
  inputIsLineIntegrals = p.Results.inputIsLineIntegrals;
  outDir = p.Results.outDir;
  printEvery = p.Results.printEvery;
  principalPoint = p.Results.principalPoint;
  projAngles = p.Results.projAngles;
  verbose = p.Results.verbose;

  if numel( dProjAngle ) == 0  &&  numel( projAngles ) == 0
    error( 'Either dProjAngle or projAngles must be provided.' );
  end

  if numel( sImg ) ~= 2  ||  min( sImg ) < 1
    error( 'sImg must be a two element array specifying the size of the output image for each slice' );
  end

  if isscalar( cropAmount ), cropAmount = [ cropAmount cropAmount ]; end

  nDetectorsY = size( sinoFlatPanel, 1 );
  nProjections = size( sinoFlatPanel, 3 );

  if numel( projAngles ) == 0
    projAngles = dProjAngle * ( 0 : nProjections - 1 );
  end

  out = [];
  if isempty( outDir )
    if isempty( cropAmount )
      out = zeros( nDetectorsY, sImg(1), sImg(2) );
    else
      out = zeros( nDetectorsY, cropAmount(1), cropAmount(2) );
    end
  end

  pbar = parforProgress( nDetectorsY );
  parfor j = 1 : nDetectorsY
    if verbose == true, pbar.progress( j, printEvery ); end   %#ok<PFBNS>

    sinoFan = squeeze( sinoFlatPanel( j, :, : ) );
    if inputIsLineIntegrals
      sinoMu = sinoFan;
    else
      if isscalar( I0 )
        sinoMu = -log( max( sinoFan, eps ) / max( I0, eps ) );
      else
        I0Slice = squeeze( I0( j, :, : ) );
        sinoMu = -log( max( sinoFan, eps ) ./ max( I0Slice, eps ) );
      end
    end
    sinoMu( ~isfinite( sinoMu ) ) = 0;

    sinoPara = fan2Parallel( sinoMu, dSizeX, f, d2Rot, 'projAngles', projAngles, ...
      'detCenter', detCenter, 'principalPoint', principalPoint );
    sinoPara( ~isfinite( sinoPara ) ) = 0;

    recon2D = filtBackProject( sinoPara, sImg, 'projAngles', projAngles, 'dSize', dSizeX );
    if numel( cropAmount ) > 0, recon2D = cropData( recon2D, cropAmount ); end

    if isempty( outDir )
      out(j,:,:) = reshape( recon2D, [ 1, sImg(:)' ] );
    else
      recon2D = max( recon2D, 0 );
      recon2D = recon2D / max( recon2D(:) );
      outFile = [ outDir, filesep, 'slice_', indx2str( j, nDetectorsY ), '.jpg' ];
      imwrite( recon2D, outFile );
    end
  end
  pbar.clean;
end
