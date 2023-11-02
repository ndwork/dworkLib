
function vol2obj( vol, varargin )
  % vol2obj( vol [, 'filename', filename, 'nTransparencies', nTransparencies, 'thresh',
  %   thresh, 'verbose', true/false );
  %
  % Written by Nicholas Dwork - Copyright 2023
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addParameter( 'filename', 'vol', @(x) true );
  p.addParameter( 'nTransparencies', 1000, @ispositive );
  p.addParameter( 'thresh', 0.05 );
  p.addParameter( 'verbose', false, @(x) islogical(x) || isnumeric(x) );
  p.parse( varargin{:} );
  filename = p.Results.filename;
  nTransparencies = p.Results.nTransparencies;
  thresh = p.Results.thresh;
  verbose = p.Results.verbose;

  if min( abs( vol(:) ) ) < 0 || max( abs( vol(:) ) ) > 1
    error( 'Assumes range of vol is [0 1]' );
  end

  objFile = [ filename, '.obj' ];   objID = fopen( objFile, 'w' );
  mtlFile = [ filename, '.mtl' ];   mtlID = fopen( mtlFile, 'w' );

  % Create the MTL file with many transparencies
  for t = 0 : nTransparencies
    vString = num2str( 1 - t / nTransparencies );
    vStringC = num2str( t / nTransparencies );
    fprintf( mtlID, [ 'newmtl trans_', num2str(t), '\n' ] );
    fprintf( mtlID, [ '  Ka ', vString, ' 0 ', vStringC, '\n' ] );
    fprintf( mtlID, [ '  Kd ', vString, ' 0 ', vStringC, '\n' ] );
    fprintf( mtlID, '  Ks 0.000 0.000 0.000\n' );
    %fprintf( mtlID, [ '  d ', num2str(t/nTransparencies), '\n' ] );
    %fprintf( mtlID, [ '  Tr ', vString, '\n\n' ] );
    fprintf( mtlID, [ '  Tr 0.9 \n\n' ] );
  end
  
  sVol = size( vol );

  fprintf( objID, '# vertices \n' );
  fprintf( objID, '\nmtllib vol.mtl\n');

  nVol = numel( vol );
  for i = 1 : nVol
    if verbose == true  &&  mod( i+1, sVol(1)*sVol(2) ) == 0
      disp([ 'Working on slice ', indx2str(z,sVol(3)), ' of ', num2str(sVol(3)) ]);
    end

    if vol( i ) < thresh, continue; end
    [y,x,z] = ind2sub( sVol, i );

    % We will make size faces for six sides of the cube
    % This will consistute 8 vertices
    fprintf( objID, [ 'v ', num2str(x),   ' ', num2str(y),   ' ', num2str(z), '\n' ] );
    fprintf( objID, [ 'v ', num2str(x+1), ' ', num2str(y),   ' ', num2str(z), '\n' ] );
    fprintf( objID, [ 'v ', num2str(x+1), ' ', num2str(y+1), ' ', num2str(z), '\n' ] );
    fprintf( objID, [ 'v ', num2str(x), ' ', num2str(y+1), ' ', num2str(z), '\n' ] );

    fprintf( objID, [ 'v ', num2str(x),   ' ', num2str(y),   ' ', num2str(z+1), '\n' ] );
    fprintf( objID, [ 'v ', num2str(x+1), ' ', num2str(y),   ' ', num2str(z+1), '\n' ] );
    fprintf( objID, [ 'v ', num2str(x+1), ' ', num2str(y+1), ' ', num2str(z+1), '\n' ] );
    fprintf( objID, [ 'v ', num2str(x),   ' ', num2str(y+1), ' ', num2str(z+1), '\n' ] );
    fprintf( objID, '\n' );

    fprintf( objID, [ 'usemtl trans_', num2str( round( vol(i) * nTransparencies ) ), '\n' ] );
    fprintf( objID, 'f -8 -7 -6 -5\n' );
    fprintf( objID, 'f -7 -3 -2 -6\n' );
    fprintf( objID, 'f -6 -2 -1 -5\n' );
    fprintf( objID, 'f -8 -7 -3 -4\n' );
    fprintf( objID, 'f -8 -4 -1 -5\n' );
    fprintf( objID, 'f -4 -3 -2 -1\n' );
    fprintf( objID, '\n' );
  end

  fclose( objID );
  fclose( mtlID );
end
