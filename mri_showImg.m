
function mri_showImg( varargin )
  % mri_showImg( [ file ] )
  %
  % Inputs:
  % file - the P*.7 file to display.  If not specified, browser window is
  %        displayed.
  %
  % Written by Nicholas Dwork

  p = inputParser;
  p.addOptional( 'file', [], @(x) true );
  p.parse( varargin{:} );
  file = p.Results.file;

  if numel( file ) == 0
    file = uigetfile('P*.7');
    if isa( file, 'double' ) && file==0, return; end
  end

  data = rawloadX( file );

  magImg = abs( fftshift( ifft2( ifftshift( data ) ) ) );

  figure;
  imshownice( magImg, 3, 'border', 0 );
end
