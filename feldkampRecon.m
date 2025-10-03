

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
  %   Either dAngle or projAngles must be supplied
  % dSizeX - the size of the detector in the dimension of rotation
  % dSizeY - the size of the detector in the dimension of the rotation axis
  % projAngles - the angle of the principal ray for each projection
  %   Either dAngle or projAngles must be supplied
  % verbose - true/false to print verbosity statements
  %
  % Written by Nicholas Dwork, Copyright 2025
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  % TODO: Density compensation may be beneficial

  p = inputParser;
  p.addParameter( 'dProjAngle', [], @isnumeric );
  p.addParameter( 'dSizeX', 1, @ispositive );
  p.addParameter( 'dSizeY', 1, @ispositive );
  p.addParameter( 'I0', 1, @ispositive );
  p.addParameter( 'printEvery', 10, @ispositive );
  p.addParameter( 'projAngles', [], @isnumeric );
  p.addParameter( 'verbose', false, @islogical );
  p.parse( varargin{:} );
  dProjAngle = p.Results.dProjAngle;
  dSizeX = p.Results.dSizeX;
  dSizeY = p.Results.dSizeY;
  I0 = p.Results.I0;
  printEvery = p.Results.printEvery;
  projAngles = p.Results.projAngles;
  verbose = p.Results.verbose;

  if numel( dProjAngle ) == 0  &&  numel( projAngles ) == 0
    error( 'Either dProjAngle or projAngles must be provided.' );
  end

  if numel( sImg ) ~= 2  ||  min( sImg ) < 1
    error( 'sImg must be a two element array specifying the size of the output image for each slice' );
  end

  nDetectorsY = size( sinoFlatPanel, 1 );
  nDetectorsX = size( sinoFlatPanel, 2 );
  nProjections = size( sinoFlatPanel, 3 );

  if numel( projAngles ) == 0
    projAngles = dProjAngle * ( 0 : nProjections-1 );
  end

  out = cell( nDetectorsY, 1, 1 );
  if verbose == true, p = parforProgress( nProjections ); end
  parfor j = 1 : nDetectorsY
    if verbose == true, p.progress( j, printEvery ); end
    sinoFan = squeeze( sinoFlatPanel(j,:,:) );
    sinoFan = -log( sinoFan / I0 );  % convert from received intensities to integrals of attenuation.
    sinoPara = fan2Parallel( sinoFan, dSizeX, f, d2Rot, 'projAngles', projAngles );
    sinoPara( ~isfinite( sinoPara ) ) = 0;
    recon2D = filtBackProject( sinoPara, sImg, 'projAngles', projAngles, 'dSize', dSizeX );
    out{j,1,1} = reshape( recon2D, [ 1, sImg(:)' ] );
  end
  if verbose == true, p.clean; end
  out = cell2mat( out );
end
