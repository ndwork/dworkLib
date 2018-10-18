
function testDworkLib
  clear; close all; rng(2);

  %% bilateralFilter
  fprintf('\nTesting bilateralFilter: \n');
  img = phantom();  sImg = size(img);
  noiseSig = 0.1;
  noise = normrnd( 0, noiseSig, sImg(1), sImg(2) );
  noisyImg = img + noise;
  denoisedImg = bilateralFilter( img, 'sigmaR', 0.1 );
  figure; imshow( [noisyImg, denoisedImg], [] );
  title('Bilateral Filter Result');
  
  %% cropData
  fprintf('\nTesting cropData: \n');
  cropped = cropData( 1:10, 5 );
  err = norm( cropped - [4 5 6 7 8], 2 );
  if err > 0, error( 'cropData test failed' ); end;
  cropped = cropData( 1:5, 4 );
  err = norm( cropped - [1 2 3 4], 2 );
  if err > 0, error( 'cropData test failed' ); end;
  cropped = cropData( 1:6, 4 );
  err = norm( cropped - [2 3 4 5] );
  if err > 0, error( 'cropData test failed' ); end;
  cropped = cropData( 1:5, 3 );
  err = norm( cropped - [2 3 4] );
  if err > 0, error( 'cropData test failed' ); end;
  disp('cropData test passed');

  %% evlautePoly
  c = [1 2 3];  x = [4 5 6];
  p1 = evaluatePoly( c, x );
  [p2,dp] = evaluatePoly( c, x );
  rightP = c(1) + c(2)*x + c(3)*x.^2;
  rightDp = c(2) + 2*c(3)*x;
  err1 = norm( p1 - rightP );
  err2 = norm( p2 - rightP );
  err3 = norm( dp - rightDp );
  if err1 + err2 + err3 < 1d-8
    disp( 'evlautePoly passed' );
  else
    error( 'evlautePoly failed' );
  end

  %% fista test 1
  M = 300;
  N = 20;
  A = rand(M,N);
  b = rand(M,1);

  bestX = A\b;

  g = @(x) 0.5 * norm( A*x - b, 2 ).^2;
  gGrad = @(x) A'*A*x - A'*b;
  proxth = @(x,t) x;  % Least squares
  x0 = zeros(N,1);
  xHat_leastSquares = fista( x0, g, gGrad, proxth );

  err = norm( xHat_leastSquares - bestX ) / M;
  if err < 1d-7
    disp( 'fista test 1 passed' );
  else
    error( 'fista test 1 failed' );
  end;

  %% fista test 2
  lambda = 50;
  M = 300;
  N = 30;
  A = rand(M,N);
  x = zeros(N,1);
  x(3)=5; x(2)=8;  x(22)=2;
  %b = A*x + rand(M,1);
  b = A*x + rand(M,1)/10;

  g = @(x) 0.5 * norm( A*x - b, 2 ).^2;
  h = @(y) norm( y, 1 );
  gGrad = @(x) A'*A*x - A'*b;
  x0 = rand(N,1);

  proxth = @(x,t) softThresh( x, lambda*t );  % Lasso
  [xHat_lasso,objValues] = fista( x0, g, gGrad, proxth, 'h', h );                          %#ok<ASGLU>
  err = norm( xHat_lasso - x, 1 ) / N;
  if err < 0.1
    disp( 'fista passed' );
  else
    error( 'fista failed' );
  end
  
  %% findDoGFeatures2D
  imgFile = '/Applications/MATLAB_R2016a.app/toolbox/images/imdata/moon.tif';
  img = double( imread( imgFile ) );
  features = findDoGFeatures2D( img, 'nFeatures', 10 );
  figure;
  showFeaturesOnImg( features, img );

  %% findRotAndTransFromPts
  nPts = 5;
  theta = pi/4;
  pts1 = 10*rand(nPts,2) - 5;
  R = [ cos(theta) -sin(theta); sin(theta) cos(theta) ];
  t = [ 5.2; -4.1 ];
  pts2 = transpose( R*transpose(pts1) ) + repmat(t', [nPts 1]);
  [R2,t2] = findRotAndTransFromPts( pts1, pts2 );
  error1 = norm( R2(:) - R(:), 2 );
  error2 = norm( t - t2, 2 );
  if error1 + error2 < 1d-7
    disp('findRotAndTransFromPts passed');
  else
    error('findRotAndTransFromPts failed');
  end

  %% findRotWithPCC
  rotation = -30;  % degrees
  img1 = padData( imread( 'cameraman.tif' ), [356 356] );
  img2 = imrotate( img1, rotation, 'crop' );
  rotation = findRotWithPCC( img1, img2 );
  rotated1 = imrotate( img1, rotation*180/pi, 'crop' );
  err = norm( img2(:) - rotated1(:) );
  if err < 1d-7
    disp('findRotWithPCC passed');
  else
    error('findRotWithPCC failed');
  end

  %% findTransRotWithPCC
  trans = [10 50];  % pixels
  rotation = -30;  % degrees
  img1 = padData( imread( 'cameraman.tif' ), [600 600] );
  tmp = circshift( img1, trans );
  img2 = imrotate( tmp, rotation, 'crop' );
  [vShift, hShift, rotation] = findTransRotWithPCC( img1, img2 );
  trans1 = imrotate( circshift( img1, [vShift hShift] ), rotation * 180/pi, ...
    'crop' );
  err = norm( img2(:) - trans1(:) );
  if err < 1d-7
    disp('findTransRotWithPCC passed');
  else
    error('findTransRotWithPCC failed');
  end

  %% findTransWithPCC
  trans = [10 200];  % pixels
  img1 = padData( imread( 'cameraman.tif' ), [356 356] );
  img2 = circshift( img1, trans );
  [vShift,hShift] = findTransWithPCC( img1, img2 );
  shifted1 = circshift( img1, [vShift hShift] );
  err = norm( img2(:) - shifted1(:) );
  if err < 1d-7
    disp('findTransWithPCC passed');
  else
    error('findTransWithPCC failed');
  end

  %% haar
  sig = rand(1,8);
  split = [1 0];
  wt = wtHaar( sig, split );
  sigHat = iwtHaar( wt, split );
  err = norm( sig - sigHat, 2 );
  if err > 1d-12
    error( ['wtHaar error: ', num2str(err,2)] );
  else
    disp('wtHaar test passed');
  end

  %% haar adjoint
  x = rand(1,8);
  y = rand(1,8);
  split = [ 1 0 ];
  wtx = wtHaar( x, split );
  iwty = iwtHaar( y, split );
  err = abs( dotP( wtx, y ) - dotP( x, iwty ) );
  if err > 1d-12
    error( 'wtHaar is not orthogonal' );
  else
    disp('wtHaar is orthogonal');
  end

  %% haar2
  im = phantom(); 
  wIm = wtHaar2( im );
  imHat = iwtHaar2( wIm );
  err = norm( im(:) - imHat(:), 2 );
  if err > 1d-12
    error(['wtHaar2 error: ', num2str(err)]);
  else
    disp('wtHaar2 passed');
  end

  %% haar2 adjoint
  x = rand(8);
  y = rand(8);
  split = [ 1 0; 0 0; ];
  wtx = wtHaar2( x, split );
  iwty = iwtHaar2( y, split );
  err = abs( dotP( wtx, y ) - dotP( x, iwty ) );
  if err > 1d-12
    error( 'wtHaar2 is not orthogonal: error' );
  else
    disp('wtHaar2 is orthogonal: passed');
  end

  %% homographyFromPts2D
  pts1 = [ [0 0]; [0 1]; [1 0]; [1 1]; ];
  H = rand(3,3);
  pts1_h = euc2Hom( pts1' );
  pts2_h = H * pts1_h;
  pts2 = hom2Euc( pts2_h )';
  dltH = homographyFromPts2D( pts1, pts2 );
  H = H ./ H(3,3);
  dltH = dltH ./ dltH(3,3);
  err = norm( H - dltH, 'fro' );
  if err < 1d-10
    disp('homographFromPts2D passed');
  else
    error('homographFromPts2D failed');
  end

  %% homographFromPts3D
  pts1 = [ [0 0 0]; [0 0 1]; [0 1 0]; [0 1 1]; [1 0 0]; [1 0 1]; ];
  H = rand(4,4);
  pts1_h = euc2Hom( pts1' );
  pts2_h = H * pts1_h;
  pts2 = hom2Euc( pts2_h )';
  dltH = homographyFromPts3D( pts1, pts2 );
  H = H ./ H(4,4);
  dltH = dltH ./ dltH(4,4);
  err = norm( H - dltH, 'fro' );
  if err < 1d-10
    disp('homographFromPts3D passed');
  else
    error('homographFromPts3D failed');
  end

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

  %% lassoCP
  M = 100;
  N = M;
  K = rand( M, N );
  x = zeros(N,1);
  x(2) = 1;
  x(3) = pi;
  x(20) = 10;
  x(N) = 8;
  b = K*x;
  gamma = 1;
  nIter = 1000;
  [xHat,residuals] = lassoCP( K, b, gamma, 'nIter', nIter );
  [xHatLS,residualsLS] = lassoCP( K, b, gamma, 'nIter', nIter );
  err = norm(xHat-xHatLS,2);
  if err > 1d-4, error( 'lassoCP failed'); end;
  disp('lassoCP passed');
  %semilogy( residuals, 'b', 'LineWidth', 2 );
  %hold on; semilogy( residualsLS, 'r', 'LineWidth', 2 );
  %legend( 'cp', 'cpls' );

  %%% lsqrFISTA
  %fprintf( '\nTesting lsqrFISTA: \n');
  %A = rand(9,8);
  %b = rand(9,1);
  %tolerance = 1d-8;
  %x = A \ b;
  %maxIter = 10000;
  %x0 = rand( 8, 1 );
  %x1 = lsqrFISTA( A, b, tolerance, maxIter, x0 );
  %err1 = norm( x1 - x, 2 ) / norm(x,2);
  %if ~isfinite(err1) || err1 > 1d-4
  %  error(['lsqrFISTA with matrix failed with error ', num2str(err1)]);
  %else
  %  disp('lsqrFISTA with matrix passed');
  %end
  %
  %function out = applyA( in, type )
  %  if nargin>1 && strcmp(type,'transp')
  %    out = A'*in;
  %  else
  %    out = A*in;
  %  end
  %end
  %x2 = lsqrFISTA( @applyA, b, tolerance, maxIter, x0 );
  %err2 = norm( x2 - x, 2 ) / norm(x,2);
  %if ~isfinite(err2) || err2 > 1d-4
  %  error(['lsqrFISTA with file handle failed with error ', num2str(err2)]);
  %else
  %  disp('lsqrFISTA with file handle passed');
  %end

  %% lsqrTikhonov
  damping = 0;
  A = rand(20,5);
  x = rand(5,1);
  b = A*x;
  xHat = lsqrTikhonov( A, b, damping );
  err = norm( x - xHat, 2 ) / norm(x,2);
  if err > 1d-10
    error(['lsqrTikhonov failed with error ', num2str(err)]);
  else
    disp('lsqrTikhonov passed');
  end

  %% lsqrTV
  gamma = 0.0;
  A = rand(20,5);
  x = rand(5,1);
  b = A*x;
  xHat = lsqrTV( A, b, gamma );
  err = norm( x - xHat, 2 ) / norm(x,2);
  if err > 1d-11
    error(['lsqrTV failed with error ', num2str(err)]);
  else
    disp('lsqrTV passed');
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

  %% matrixVolProd
  A = rand(3,4);
  vol = rand(4,5,10000);
  out1 = zeros( [ size(A,1) size(vol,2) size(vol,3) ] );
  for m=1:size(vol,3)
    out1(:,:,m) = A * vol(:,:,m);
  end
  out2 = matrixVolProd( A, vol );
  err = max( abs( out1(:) - out2(:) ) );
  if err > 1d-14
    error(['matrixVolProd failed with error ', num2str(err)]);
  else
    disp('matrixVolProd passed');
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

  %% padData
  fprintf('\nTesting padData: \n');
  padded = padData( [1 2 3], 5 );
  err = norm( padded - [0 1 2 3 0] );
  if err>0, error(['padData failed with error ', num2str(err)]); end;
  padded = padData( [1 2 3], 4 );
  err = norm( padded - [0 1 2 3] );
  if err>0, error(['padData failed with error ', num2str(err)]); end;
  padded = padData( [1 2], 4 );
  err = norm( padded - [0 1 2 0] );
  if err>0, error(['padData failed with error ', num2str(err)]); end;
  padded = padData( [1 2], 3 );
  err = norm( padded - [1 2 0] );
  if err>0, error(['padData failed with error ', num2str(err)]); end;

  % test is cropData and padData are adjoints of each other
  x = rand(5,1);  y = rand(3,1);
  err = dotP( cropData(x,3), y ) - dotP( x, padData(y,5) );
  if err>0, error(['padData failed with error ', num2str(err)]); end;
  x = rand(6,1); y = rand(4,1);
  err = dotP( cropData(x,4), y ) - dotP( x, padData(y,6) );
  if err>0, error(['padData failed with error ', num2str(err)]); end;
  x = rand(5,1); y = rand(4,1);
  err = dotP( cropData(x,4), y ) - dotP( x, padData(y,5) );
  if err>0, error(['padData failed with error ', num2str(err)]); end;
  x = rand(6,1); y = rand(5,1);
  err = dotP( cropData(x,5), y ) - dotP( x, padData(y,6) );
  if err>0, error(['padData failed with error ', num2str(err)]); end;
  disp('padData passed');

  %% parforProgress
  N = 10;
  p = parforProgress( N );
  parfor n=1:N
    p.progress( n );
    pause(rand); % Replace with real code
  end
  p.clean;

  %% plotnice
  x = 1:10; y = 2*x;
  figure; subplot(4,1,1); plotnice(x);
  subplot(4,1,2); plotnice(x,y);
  subplot(4,1,3); plotnice(x, y, 'ro');
  subplot(4,1,4); plotnice( x, y, 'g+', 'MarkerSize', 20 );
  disp('Plotnice passed');

  %% powerIteration
  fprintf('\nTesting powerIteration: \n');
  M = rand(3);
  est = powerIteration( M, false );
  normM = norm( M );
  err = abs( est - normM ) / normM;
  if err > 1d-5, error('Power Iteration failed'); end;

  symmM = M + M';
  [est,flag] = powerIteration( symmM, true );                                              %#ok<ASGLU>
  normSymmM = norm( symmM );
  err = abs( est - normSymmM ) / normM;
  if err > 1d-5, error('Power Iteration failed'); end;

  function out = applyM( x, type )
    if (nargin>1) && strcmp( type, 'transp')
      out = M'*x;
    else
      out = M*x;
    end
  end
  x0 = rand( size(M,2), 1 );
  est = powerIteration( @applyM, false, x0 );
  err = abs( est - normM ) / normM;
  if err > 1d-5, error('Power Iteration failed'); end;
  disp('powerIteration passed');

  function out = applyMsymm( x )
    out = symmM * x;
  end
  x0 = rand( size(symmM,2), 1 );
  est = powerIteration( @applyMsymm, true, x0 );
  err = abs( est - normSymmM ) / normSymmM;
  if err > 1d-5, error('Power Iteration failed'); end;

  disp('powerIteration passed');

  %% ransacDltHomographyFromPts2D
  pts1 = [ [0 0]; [0 1]; [1 0]; [1 1]; ];
  pts1 = [ pts1; rand(5,2); ];
  H = rand(3,3);
  pts1_h = euc2Hom( pts1' );
  pts2_h = H * pts1_h;
  pts2 = hom2Euc( pts2_h )';
  dltH = ransacDltHomographyFromPts2D( pts1, pts2, 1 );
  H = H ./ H(3,3);
  dltH = dltH ./ dltH(3,3);
  err = norm( H - dltH, 'fro' );
  if err < 1d-8
    disp('ransacDltHomographyFromPts2D passed');
  else
    error('ransacDltHomographyFromPts2D failed');
  end

  %% ransacDltHomographyFromPts3D
  pts1 = [ [0 0 0]; [0 0 1]; [0 1 0]; [0 1 1]; [1 0 0]; [1 0 1]; ];
  pts1 = [ pts1; rand(5,3); ];
  H = rand(4,4);
  pts1_h = euc2Hom( pts1' );
  pts2_h = H * pts1_h;
  pts2 = hom2Euc( pts2_h )';
  dltH = ransacDltHomographyFromPts3D( pts1, pts2, 1 );
  H = H ./ H(4,4);
  dltH = dltH ./ dltH(4,4);
  err = norm( H - dltH, 'fro' );
  if err < 1d-8
    disp('ransacDltHomographyFromPts3D passed');
  else
    error('ransacDltHomographyFromPts3D failed');
  end

  %% ransacRotAndTrans
  nPts = 5;
  theta = pi/4;
  distThresh = 5;
  pts1 = 10*rand(nPts,2) - 5;
  R = [ cos(theta) -sin(theta); sin(theta) cos(theta) ];
  t = [ 5.2; -4.1 ];
  pts2 = transpose( R*transpose(pts1) ) + repmat(t', [nPts 1]);
  [R2,t2] = ransacRotAndTrans( pts1, pts2, distThresh );
  error1 = norm( R2(:) - R(:), 2 );
  error2 = norm( t - t2, 2 );
  if error1 + error2 < 1d-7
    disp('ransacRotAndTrans passed');
  else
    error('ransacRotAndTrans failed');
  end

  %% rootsOfQuadratic
  a=rand(1) * sign( rand(1) - 0.5 );
  b=rand(1) * sign( rand(1) - 0.5 );
  c=rand(1) * sign( rand(1) - 0.5 );
  [x1,x2] = rootsOfQuadratic( a, b, c );
  v1 = a*x1*x1 + b*x1 + c;
  v2 = a*x2*x2 + b*x2 + c;
  if( abs(v1) < 1d-7 && abs(v2) < 1d-7 )
    disp('rootsOfQuadratic passed');
  else
    error('rootsOfQuadratic failed');
  end

  %% shearImg
  imgFile = '/Applications/MATLAB_R2016a.app/toolbox/images/imdata/moon.tif';
  img = double( imread( imgFile ) );
  ySheared = shearImg( img, pi/8 );
  xSheared = shearImg( img, pi/8, 2 );
  figure; imshow( [ img xSheared ySheared ], [] );
  title('original / xSheared / ySheared');
  

  %% showLibs
  fprintf( '\nTesting showLibs: \n');
  showLibs();

  %% showLibFiles
  fprintf('\nTesting showLibFiles: \n');
  showLibFiles( 'dworkLib' );

  %% vecVolMatrixProd
  v = rand(1,5);
  vol = rand(5,3,100000);
  sVol = size( vol );
  out1 = zeros( 1, sVol(2), sVol(3) );
  for m = 1:sVol(3)
    out1(:,:,m) = v * vol(:,:,m);
  end
  out2 = vecVolMatrixProd( v, vol );
  err = max( abs( out1(:) - out2(:) ) );
  if err < 1d-10
    disp('vecVolMatrixProd passed');
  else
    error('vecVolMatrixProd failed');
  end

  %% volMatrixProd
  vol = rand(3,4,10000);
  A = rand(4,2);
  out1 = zeros( [ size(vol,1) size(A,2) size(vol,3) ] );
  for m=1:size(vol,3)
    out1(:,:,m) = vol(:,:,m) * A;
  end
  out2 = volMatrixProd( vol, A );

  err = max( abs( out1(:) - out2(:) ) );
  if err < 1d-10
    disp('volMatrixProd passed');
  else
    error('volMatrixProd failed');
  end


  %% volVecMatrixProd
  vol = rand(3,4,100000);
  v = rand(4,1);

  sVol = size( vol );
  out1 = zeros( sVol(1), 1, sVol(3) );
  for m = 1:sVol(3)
    out1(:,m) = vol(:,:,m) * v;
  end
  out2 = volVecMatrixProd( vol, v );
  err = max( abs( out1(:) - out2(:) ) );
  if err < 1d-7
    disp('volVectorMatrixProd passed');
  else
    error('volVectorMatrixProd failed');
  end


  %% volVecsMatrixProd
  vol = rand(3,4,100000);
  vs = rand(4,100000);

  sVol = size( vol );
  out1 = zeros( sVol(1), 1, sVol(3) );
  for m = 1:sVol(3)
    out1(:,m) = vol(:,:,m) * vs(:,m);
  end
  out2 = volVecsMatrixProd( vol, vs );
  err = max( abs( out1(:) - out2(:) ) );
  if err < 1d-7
    disp('volVectorMatrixProd passed');
  else
    error('volVectorMatrixProd failed');
  end


  %% volVolMatrixProd
  M = rand(3,3);
  vol1 = rand(3,3,100000);
  vol2 = rand( size(vol1) );
  tmp = zeros(size(vol1));
  for m=1:size(vol1,3)
    tmp(:,:,m) = vol1(:,:,m) * vol2(:,:,m);
  end
  asdf = volVolMatrixProd( vol1, vol2 );
  err = max( abs( asdf(:) - tmp(:) ) );
  if err < 1d-7
    disp('volVolMatrixProd passed');
  else
    error('volVolMatrixProd failed');
  end

end
