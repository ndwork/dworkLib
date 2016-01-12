
function testDworkLib
  clear; close all; rng(1);

  %% bilateralFilter
  fprintf('\nTesting bilateralFilter: \n');
  img = phantom();  sImg = size(img);
  noiseSig = 0.1;
  noise = normrnd( 0, noiseSig, sImg(1), sImg(2) );
  noisyImg = img + noise;
  denoisedImg = bilateralFilter( img, 'sigmaR', 0.1 );
  figure; imshow( [noisyImg, denoisedImg], [] );
  title('Bilateral Filter Result');

  %% isEven - 1D data
  fprintf( '\nTesting isEven (1D): \n');
  A1 = rand(5,1);
  A1odd = 0.5 * ( A1 + flipud(A1) );
  if ~isEven( A1odd )
    error('isEven (1D) failed');
  else
    disp('isEven (1D) passed');
  end
  
  %% isEven - 2D data
  fprintf( '\nTesting isEven (2D): \n');
  A1 = rand(5,5);
  A1odd = 0.5 * ( A1 + rot90(A1,2) );
  if ~isEven( A1odd )
    error('isEven (2D) failed');
  else
    disp('isEven (2D) passed');
  end
  
  %% isHermitian - 1D data
  fprintf( '\nTesting isHermitian (1D): \n');
  A1 = rand(5,1);
  A1even = 0.5 * ( A1 + flipud(A1) );
  A1odd = 0.5 * ( A1 - flipud(A1) );
  A1hermitian = A1even + 1i * A1odd;
  if ~isHermitian( A1hermitian )
    error('isEven (1D) failed');
  else
    disp('isEven (1D) passed');
  end

  %% isHermitian - 2D data
  fprintf( '\nTesting isHermitian (2D): \n');
  A1 = rand(5,1);
  A1even = 0.5 * ( A1 + rot90(A1,2) );
  A1odd = 0.5 * ( A1 - rot90(A1,2) );
  A1hermitian = A1even + 1i * A1odd;
  if ~isHermitian( A1hermitian )
    error('isEven (2D) failed');
  else
    disp('isEven (2D) passed');
  end

  %% isOdd - 1D data
  fprintf( '\nTesting isOdd (1D): \n');
  A1 = rand(5,1);
  A1odd = 0.5 * ( A1 - flipud(A1) );
  if ~isOdd( A1odd )
    error('isOdd (1D) failed');
  else
    disp('isOdd (1D) passed');
  end

  %% isOdd - 2D data
  fprintf( '\nTesting isOdd (2D): \n');
  A1 = rand(5,5);
  A1odd = 0.5 * ( A1 - rot90(A1,2) );
  if ~isOdd( A1odd )
    error('isOdd (2D) failed');
  else
    disp('isOdd (2D) passed');
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

  %% nonlocal mean
  fprintf('\nTesting nonlocal means: \n');
  img = phantom();  sImg = size(img);
  noiseSig = 0.1;
  noise = normrnd( 0, noiseSig, sImg(1), sImg(2) );
  noisyImg = img + noise;
  denoisedImg = nonlocalMeans( img, 'sigmaS', 0.02 );
  figure; imshow( [noisyImg, denoisedImg], [] );
  title('Nonlocal Means Result');
  
  %% showLibs
  fprintf( '\nTesting showLibs: \n');
  showLibs();

  %% showLibFiles
  fprintf('\nTesting showLibFiles: \n');
  showLibFiles( 'dworkLib' );

end
