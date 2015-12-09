
function testDworkLib
  clear; close all; rng(1);

  %% isEven
  fprintf( '\nTesting isEven: \n');
  A1 = rand(5,5);
  A1even = 0.5 * ( A1 + rot90(A1,2) );
  if ~isEven( A1even )
    error('isEven failed');
  else
    disp('isEven passed');
  end

  %% makeDftMatrix
  fprintf( '\nTesting makeDftMatrix: \n');
  M = 100;
  standard = fft( eye(M) );
  out = makeDftMatrix( M, M );
  diff = out - standard;
  err = max( abs( diff(:) ) );
  if err > 1d-14
    error(['makeDftMatrix failed with error ', num2str(err)]);
  else
    disp('makeDftMatrix passed');
  end

  %% showLibs
  fprintf( '\nTesting showLibs: \n');
  showLibs();

  %% showLibFiles
  fprintf('\nTesting showLibFiles: \n');
  showLibFiles( 'dworkLib' );

end
