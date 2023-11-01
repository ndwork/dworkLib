
function vol2obj( vol, varargin )

  p = inputParser;
  p.addParameter( 'filename', 'vol', @(x) true );
  p.parse( varargin{:} );
  filename = p.Results.filename;

  if min( abs( vol(:) ) ) < 0 || max( abs( vol(:) ) ) > 1
    error( 'Assumes range of vol is [0 1]' );
  end

  objFile = [ filename, '.obj' ];   objID = fopen( objFile, 'w' );
  mtlFile = [ filename, '.mtl' ];   mtlID = fopen( mtlFile, 'w' );

  % Create the MTL file with many transparencies
  nTransparencies = 10;
  for t = 0 : nTransparencies
    fprintf( mtlID, [ 'newmtl trans_', num2str(t), '\n' ] );
    fprintf( mtlID, '  Ka 1.000 1.000 1.000\n' );
    fprintf( mtlID, '  Kd 1.000 1.000 1.000\n' );
    fprintf( mtlID, '  Ks 0.000 0.000 0.000\n' );
    fprintf( mtlID, [ '  d ', num2str(1-t/nTransparencies), '\n' ] );
    fprintf( mtlID, [ '  Tr ', num2str(t/nTransparencies), '\n\n' ] );
  end

  sVol = size( vol );  
  faces = zeros( [ sVol 6 4 ] );

  fprintf( objID, '# vertices \n' );

  nVol = numel( vol );
nVol = 2
  for i = 1 : nVol
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

    % Now we will store an array of faces for later
    v = (i-1) * 8 + 1;
    faces(y,x,z,1,:) = [ v   v+1 v+2 v+3 ];
    faces(y,x,z,2,:) = [ v+1 v+5 v+6 v+2 ];
    faces(y,x,z,3,:) = [ v+2 v+6 v+7 v+3 ];
    faces(y,x,z,4,:) = [ v   v+1 v+5 v+4 ];
    faces(y,x,z,5,:) = [ v   v+4 v+7 v+3 ];
    faces(y,x,z,6,:) = [ v+4 v+5 v+6 v+7 ];
  end

  fprintf( objID, '\n' );
  fprintf( objID, '# faces \n' );
  for i = 1 : nVol
    [y,x,z] = ind2sub( sVol, i );
    fprintf( objID, [ 'usemtl trans_', num2str( round( vol(i) * nTransparencies ) ), '\n' ] );
    for face = 1 : 6
      faceString = arrayfun( @num2str, faces(y,x,z,face,:), 'UniformOutput', 0 );
      fprintf( objID, [ 'f ', strjoin( faceString, ', ' ), '\n' ] );
    end
    fprintf( objID, '\n' );
  end

  fclose( objID );
  fclose( mtlID );
end
