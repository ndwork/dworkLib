
function testDworkLib
  clear; close all; rng(2);

  %% admm
  fprintf( '\nTesting admm: \n' );
  %minimize ||Ax - b||_2^2 + lambda ||x||+1fa
  A = rand( 7, 3 ) * 100;
  x = rand( 3, 1 ) * 100;
  x0 = zeros( size( x ) );
  b = A * x;
  lambda = 15;
  f = @(x) lambda * norm( x, 1 );
  g = @(x) 0.5 * norm( x - b, 2 )^2;
  proxf = @(x,t) softThresh( x, lambda * t );
  proxg = @(y,t) proxL2Sq( y, 1, b );  % TODO:  should there be a t input here?
  t = 1d-3;  % admm step size
  [xStar,objValues] = admm( x0, proxf, proxg, t, 'A', A, 'f', f, 'g', g, 'N', 1000 );   %#ok<ASGLU>
  fprintf('\nadmm ran to completion. \n');

  %% applyC_2D and applyCT_2D
  fprintf( '\nTesting applyC_2D and applyCT_2D are adjoints: \n' );
  nSpokes = 60;
  nPtsPerSpoke = 30;
  traj = mri_makeTrajPts( 2, 'radial', nSpokes, nPtsPerSpoke );

  nD = 2;
  N = [ 128 128 ];

  Ny = 2 * N(1);   [kCy,Cy,~] = makeKbKernel( Ny, Ny );
  Nx = 2 * N(2);   [kCx,Cx,~] = makeKbKernel( Nx, Nx );

  C_traj_N = @(x) applyC_2D( x, traj, N, kCy, kCx, Cy, Cx );
  CT_traj_N = @(x) applyCT_2D( x, traj, N, kCy, kCx, Cy, Cx );

  F = rand( size(traj,1), nD ) + 1i * rand( size( traj, 1 ), nD );
  [ checkC_traj_N, errCheckC ] = checkAdjoint( F, C_traj_N, 'fAdj', CT_traj_N );
  if checkC_traj_N == false
    error([ 'applyC / applyCT failed with error: ', num2str( errCheckC ) ]);
  end

  C_N_traj = @(x) applyC_2D( x, N, traj, kCy, kCx, Cy, Cx );
  CT_N_traj = @(x) applyCT_2D( x, N, traj, kCy, kCx, Cy, Cx );
  
  F = rand( [ N(1) N(2) nD ] ) + 1i * rand( [ N(1) N(2) nD ] );
  [ checkC_N_traj, errCheckC ] = checkAdjoint( F, C_N_traj, 'fAdj', CT_N_traj );
  if checkC_N_traj == false
    error([ 'applyC / applyCT failed with error: ', num2str( errCheckC ) ]);
  end

  nReadout = 50;
  nLines = 9;
  dkLine = 0.02;
  nAngles = 45;
  newTraj = mri_makeTrajPts( 2, 'propeller', nReadout, nLines, dkLine, nAngles );

  C_traj_newTraj = @(x) applyC_2D( x, traj, newTraj, kCy, kCx, Cy, Cx );
  CT_traj_newTraj = @(x) applyCT_2D( x, traj, newTraj, kCy, kCx, Cy, Cx );

  F = rand( size( traj, 1 ), nD ) + 1i * rand( size( traj, 1 ), nD );  
  [ checkC_N_traj, errCheckC ] = checkAdjoint( F, C_traj_newTraj, 'fAdj', CT_traj_newTraj );
  if checkC_N_traj == false
    error([ 'applyC failed with error: ', num2str( errCheckC ) ]);
  end
  
  disp( 'applyC and applyCT passed' );


  %% bilateralFilter
  fprintf('\nTesting bilateralFilter: \n');
  img = phantom();  sImg = size(img);
  noiseSig = 0.1;
  noise = normrnd( 0, noiseSig, sImg(1), sImg(2) );
  noisyImg = img + noise;
  denoisedImg = bilateralFilter( img, 'sigmaR', 0.1 );
  figure; imshowscale( [noisyImg, denoisedImg], 3 );
  title('Bilateral Filter Result');

  %% bilinInterp2
  rng(6);
  x = rand( 4, 1 );  x = sort( x );
  y = rand( 3, 1 );  y = sort( y );
  v = rand( numel(y), numel(x) );
  xq = sort( rand(5,1) );  xq = [ xq; rand(1);  max(x); max(x); ];
  yq = sort( rand(5,1) );  yq = [ yq;  max(y); rand(1); max(y); ];
  u = rand( numel(yq), 1 );

  Av = bilinInterp2( x, y, v, xq, yq );
  Av_correct = interp2( x, y, v, xq, yq, 'bilinear', 0 );
  errInterp = norm( Av(:) - Av_correct(:) ) / norm( Av(:) );
  if errInterp > 1d-14
    error([ 'bilinInterp2 failed with error ', num2str( errInterp ) ]);
  end
  ATu = bilinInterp2( x, y, u, xq, yq, 'op', 'transp' );
  dp1 = dotP( Av,   u );
  dp2 = dotP(  v, ATu );
  err = abs( dp1 - dp2 ) / abs( dp1 );
  if err > 1d-14
    error([ 'bilinInterp2 adjoint failed with error ', num2str( errInterp ) ]);
  else
    disp( 'bilinInterp2 passed' );
  end

  %% binarySearch
  fprintf( '\nTesting binarySearch: \n' );
  trueRoot = 3;
  f = @(x) 2 * ( x - trueRoot );
  xLB = -10;  xUB = 10;
  tol = 1d-4;
  xRoot = binarySearch( f, xLB, xUB, 'tol', tol, 'nMax', 10000 );
  if abs( xRoot - trueRoot ) < tol
    disp( 'binarySearch passed' );
  else
    error( 'binarySearch failed' );
  end

  %% bisection method
  fprintf('\nTesting bisection method: \n');
  myFunc = @(x) x^2 - 80;
  myFuncRoot = bisectionMethod( myFunc, 0, 100, 'nMax', 10000 );
  if myFunc( myFuncRoot ) < 1d-12
    disp( 'bisection passed' );
  else
    error( 'bisection failed' );
  end

  %% checkProxConj
  x = rand(10,1) + 1i * rand(10,1);
  t = 0.2;
  sigma = 0.3;
  myProx = @(y,t) proxL1Complex( y, t );
  myProxConj = @(y,sigma,t) proxConjL1( y, sigma, t );
  [proxCheck,proxErr] = checkProxConj( x, myProx, myProxConj, 'sigma', sigma, 't', t );
  if proxCheck == true
    disp( 'checkProxConj passed' );
  else
    error([ 'checkProxConj failed with err: ', num2str(proxErr) ]);
  end

  %% circConv
  a = rand( 64, 128 );
  k = rand( 6, 7 );
  out1 = cropData( conv2( a, k, 'same' ), size(k) );
  out2 = cropData( circConv( k, a ), size(k) );
  circConvErr = norm( out1(:) - out2(:) ) / norm( out1(:) );
  if circConvErr < 1d-8
    disp( 'circConv test passed' );
  else
    error([ 'circConv test failed with err: ', num2str(circConvErr) ]);
  end

  %% circConv adjoint
  img = rand( 64, 128 );
  filt = rand( 6, 7 );

  blur = @(x) circConv( x, filt );
  blurT = @(x) circConv( x, filt, 'transp' );

  [blurCheck,checkErr] = checkAdjoint( img, blur, 'fAdj', blurT );
  if blurCheck == false
    error([ 'circConv adjoint test failed with error: ', num2str( checkErr ) ]);
  end
  disp( 'circConv test passed' );

  %% contoursToPolyhedron
  contours = cell(3, 1);
  contours{1} = [[0 0]; [0 1]; [1 1];  [0.7 0.7]; [1 0];];
  contours{2} = [[0 0.5]; [0.1 1.5]; [1 1.4]; [1.1 0.5];];
  contours{3} = [[0 0]; [0.25 0.25]; [0 1]; [1 1]; [1 0];];
  contours{4} = [[0 0.5]; [0.1 1.5]; [1 1.4]; [1.1 0.5];];
  triangles = contoursToPolyhedron( contours );
  plotTriangles( triangles );

  %% cropData
  fprintf('\nTesting cropData: \n');
  cropped = cropData( 1:10, 5 );
  err = norm( cropped - [4 5 6 7 8], 2 );
  if err > 0, error( 'cropData test failed' ); end
  cropped = cropData( 1:5, 4 );
  err = norm( cropped - [1 2 3 4], 2 );
  if err > 0, error( 'cropData test failed' ); end
  cropped = cropData( 1:6, 4 );
  err = norm( cropped - [2 3 4 5] );
  if err > 0, error( 'cropData test failed' ); end
  cropped = cropData( 1:5, 3 );
  err = norm( cropped - [2 3 4] );
  if err > 0, error( 'cropData test failed' ); end
  disp('cropData test passed');

  %% crossCorrelate
  img = rand( 64, 128 );
  filt = rand( 6, 7 );

  blur = @(x) crossCorrelate( x, filt );
  blurT = @(x) crossCorrelate( x, filt, 'transp' );

  [blurCheck,checkErr] = checkAdjoint( img, blur, 'fAdj', blurT );
  if blurCheck == false
    error([ 'crossCorrelate test failed with error: ', num2str( checkErr ) ]);
  end
  disp( 'crossCorrelate test passed' );
  
  %% deaubechies
  sig = rand(8,1);
  split = [1 0];
  wt = wtDaubechies( sig, split );
  sigHat = iwtDaubechies( wt, split );
  err = norm( sig - sigHat, 2 );
  if err > 1d-12
    error( ['wtDaubechies error: ', num2str(err,2)] );
  else
    disp('wtDaubechies test passed');
  end

  %% deaubechies adjoint
  x = rand(8,1);
  y = rand(8,1);
  split = [ 1 0 ];
  wtx = wtDaubechies( x, split );
  iwty = iwtDaubechies( y, split );
  err = abs( dotP( wtx, y ) - dotP( x, iwty ) );
  if err > 1d-12
    error( 'wtHaar is not orthogonal' );
  else
    disp('wtHaar looks to be orthogonal');
  end

  %% dworkLib
  dworkLib

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

  %% evlautePoly2
  x = 4;
  y = 5;
  c = reshape( 1:4, [2 2] );
  p = evaluatePoly2( c, x, y );
  rightP = 0;
  for i = 1:size(c,2)
    for j = 1:size(c,1)
      rightP = rightP + c(j,i) * x.^(i-1) * y.^(j-1);
    end
  end
  err = norm( p - rightP );
  if err < 1d-8
    disp( 'evlautePoly2 passed' );
  else
    error([ 'evlautePoly2 failed with error, ', num2str(err) ]);
  end

  %% fftnh
  x = rand( 5, 7, 6, 9 );
  [fftnCheckOut,fftnCheckErr] = checkAdjoint( x, @fftn, @fftnh );
  if fftnCheckOut == true
    disp( 'fftnh check passed' );
  else
    disp([ 'fftn check failed with error: ', num2str(fftnCheckErr) ]);
  end

  %% findDoGFeatures2D
  imgFile = '/Applications/MATLAB_R2019b.app/toolbox/images/imdata/moon.tif';
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

  %% fista test 1
  M = 300;
  N = 20;
  A = rand(M,N);
  b = rand(M,1);

  bestX = A \ b;

  %g = @(x) 0.5 * norm( A*x - b, 2 ).^2;
  gGrad = @(x) A' * A * x - A' * b;
  proxth = @(x,t) x;  % Least squares
  x0 = zeros(N,1);
  
  normATA = norm( A' * A );
  t = 0.999 / normATA;
  xHat_leastSquares = fista( x0, gGrad, proxth, 't', t, 'tol', 1d-5, 'N', 1000 );
  
  err = norm( xHat_leastSquares - bestX ) / M;
  if err < 1d-6
    disp( 'fista test 1 passed' );
  else
    error( 'fista test 1 failed' );
  end

  %% fista test 2
  lambda = 50;
  M = 300;
  N = 30;
  A = rand(M,N);
  x = zeros(N,1);
  x(3)=5; x(2)=8;  x(22)=2;
  %b = A*x + rand(M,1);
  b = A*x + rand(M,1) / 10;

  g = @(x) 0.5 * norm( A*x - b, 2 ).^2;
  h = @(y) norm( y, 1 );
  gGrad = @(x) A'*A*x - A'*b;
  x0 = rand(N,1);

  proxth = @(x,t) softThresh( x, lambda * t );  % Lasso

  normATA = norm( A' * A );
  t = 0.999 / normATA;

  [ xHat_lasso, objValues ] = fista( x0, gGrad, proxth, 't', t, 'g', g, 'h', h, 'N', 10000 );  %#ok<ASGLU>

  err = norm( xHat_lasso - x, 1 ) / N;
  if err < 0.1
    disp( 'fista test 2 passed' );
  else
    error( 'fista test 2 failed' );
  end

  
  %% fista_wLS  ( with line search )
  M = 300;
  N = 20;
  A = rand(M,N);
  b = rand(M,1);

  bestX = A \ b;

  g = @(x) 0.5 * norm( A*x - b, 2 ).^2;
  gGrad = @(x) A' * A * x - A' * b;
  proxth = @(x,t) x;  % Least squares
  x0 = zeros(N,1);

  %normATA = norm( A' * A );
  %t = 0.999 / normATA;
  t = 1;
  xHat_leastSquares = fista_wLS( x0, g, gGrad, proxth, 't0', t, 'tol', 1d-6, 'N', 1000, 'verbose', true );

  err = norm( xHat_leastSquares - bestX ) / M;
  if err < 1d-7
    disp( 'fista test 1 passed' );
  else
    error( 'fista test 1 failed' );
  end

  %% fista_wRestart
  M = 300;
  N = 20;
  A = rand(M,N);
  b = rand(M,1);

  bestX = A \ b;

  %g = @(x) 0.5 * norm( A*x - b, 2 ).^2;
  gGrad = @(x) A' * A * x - A' * b;
  proxth = @(x,t) x;  % Least squares
  x0 = zeros(N,1);

  normATA = norm( A' * A );
  t = 0.999 / normATA;
  xHat_leastSquares = fista_wRestart( x0, gGrad, proxth, 't', t, 'tol', 1d-6, 'N', 1000 );

  err = norm( xHat_leastSquares - bestX ) / M;
  if err < 1d-7
    disp( 'fista_wRestart passed' );
  else
    error( 'fista_wRestart failed' );
  end

  %% fitPolyToData2
  xOrder = 3;
  yOrder = 2;
  nPtsX = 5;
  nPtsY = 3;

  c = rand( yOrder+1, xOrder+1 ) * 10;

  x = linspace( 1, 10, nPtsX );
  y = linspace( 1, 10, nPtsY );
  [x,y] = meshgrid( x, y );
  x = x(:);  y = y(:);

  z = evaluatePoly2( c, x, y );
  cFit = fitPolyToData2( xOrder, yOrder, x(:), y(:), z );

  err = norm( c(:) - cFit(:) );
  if err > 1d-8
    error([ 'fitPolyToData2 failed with error ', num2str(err) ]);
  else
    disp( 'fitPolyToData2 passed' );
  end

  %% flipAboutIndx
  in = 1:5;
  flipped = flipAboutIndx( in, [1 2] );
  err = norm( flipped - [ 3 2 1 5 4 ] );
  if err > 0
    error([ 'flipAboutIndx failed with error ', num2str(err) ]);
  else
    disp( 'flipAboutIndx passed' );
  end
  
  %% goldenSectionSearch
  fprintf( '\nTesting goldenSectionSearch: \n' );
  trueMin = 0.7809;
  f = @(x) x^4 - 14 * x^3 + 60 * x^2 - 70 * x;
  xLB = 0;  xUB = 2;
  tol = 1d-4;
  xMin = goldenSectionSearch( f, xLB, xUB, 'tol', tol, 'nMax', 10000 );
  if abs( xMin - trueMin ) < tol
    disp( 'goldenSectionSearch passed' );
  else
    error( 'goldenSectionSearch failed' );
  end

  %% gradDescent
  fprintf( '\nTesting gradDescent: \n' );
  A = rand( 10, 3 );
  x = rand( 3, 1 );
  b = A * x;

  g = @(in) 0.5 * norm( A * in - b ).^2;
  gGrad = @(in) A' * A * in - A' * b;

  x0 = zeros( size( x ) );
  [xStar,oValues,relDiffs] = gradDescent( x0, gGrad, 'g', g, 't', 1d-2, 'N', 10000 );

  err = norm( xStar - x );
  if err < 1d-10
    disp( 'gradDescent passed' );
  else
    disp( 'gradDescent failed' );
  end


  %% haar
  sig = rand(8,1);
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
  x = rand(8,1);
  y = rand(8,1);
  split = [ 1 0 ];
  wtx = wtHaar( x, split );
  iwty = iwtHaar( y, split );
  err = abs( dotP( wtx, y ) - dotP( x, iwty ) );
  if err > 1d-12
    error( 'wtHaar is not orthogonal' );
  else
    disp('wtHaar looks to be orthogonal');
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


  %% iGrid_2D - 2D data
  N = [ 256 256 ];
  p = 0.1;

  % Let the image be sinc( x / W )
  imgCoords = size2imgCoordinates( N );
  [xs, ys] = meshgrid( imgCoords{2}, imgCoords{1} );
  fx = p * sinc( p * xs );
  fy = p * sinc( p * ys );
  f = fx .* fy;
  figure; imshowscale( f, 5 );

  % The inverse Fourier transform is W * rect( W x )
  % Make the right answer for the Fourier data
  nK = 20000;
  kTraj = rand( nK, 2 ) - 0.5;
  % Fourier values of a rect function with width p
  F = ( abs( kTraj(:,1) ) < ( 0.5 * p ) )  &  ...
      ( abs( kTraj(:,2) ) < ( 0.5 * p ) );
  F = F * 1.0;

  figure;  scatterIntensity( kTraj(:,1), kTraj(:,2), F );
  axis( [ -0.5 0.5 -0.5 0.5 ] );    axis equal;
  colorbarnice;  titlenice( 'Right answer' );

  % Determine the estimate and compare
  F_hat = iGrid_2D( f, kTraj );
  figure;  scatterIntensity( kTraj(:,1), kTraj(:,2), abs(F_hat) );
  axis( [ -0.5 0.5 -0.5 0.5 ] );  axis equal;
  colorbarnice;  titlenice( 'Estimate' );

  err = abs( F - F_hat );
  figure;  scatterIntensity( kTraj(:,1), kTraj(:,2), abs( err ) );
  axis( [ -0.5 0.5 -0.5 0.5 ] );  axis equal;
  colorbarnice;  titlenice( '|err|' );

  relErr = norm( F(:) - F_hat(:) ) / norm( F(:) );
  disp( relErr );
  if relErr >= 0.2
    error([ 'iGrid_2D failed with relative err: ', num2str(relErr) ]);
  else
    disp( 'iGrid_2D passed' );
  end


  %% Make sure iGrid_2D and iGridT_2D are adjoints
  sizeX = [ 256 256 ];
  nY = 50000;
  nD = 1;
  kTraj = rand( nY, 2 ) - 0.5;
  x = rand( [ sizeX nD ] );
  y = rand( nY, nD );
  Ax = iGrid_2D( x, kTraj );
  dp1 = dotP( Ax, y );

  ATy = iGridT_2D( y, kTraj, sizeX );
  dp2 = dotP( x, ATy );
  err = abs( dp1 - dp2 ) / abs( dp1 );
  disp([ 'iGrid/iGridT 2D Adjointness error:  ', num2str(err) ]);

  
  %% Make sure iGrid_2D and grid_2D approximately undo each other
  nSpokes = 360;
  nPtsPerSpoke = 150;
  sImg = [ 128 128 ];
  traj = mri_makeTrajPts( 2, 'radial', nSpokes, nPtsPerSpoke );
  xRect = rect( 1:sImg(2), 25 );
  yRect = rect( 1:sImg(1), 25 );
  img = yRect(:) * xRect(:)';
  img = circshift( img, [ 20 10 ] );
  FVals = iGrid_2D( img, traj );
  weights = makePrecompWeights_2D( traj, 'sImg', size( img ) );
  gridded = grid_2D( FVals, traj, size( img ), weights );
  relErr = norm( img(:) - gridded(:) ) / norm( img(:) );
  if relErr >= 0.05
    error([ 'iGrid_2D and grid_2D do not approximately undo each other with relative err: ', ...
      num2str(relErr) ]);
  else
    disp( 'iGrid_2D / grid_2D passed' );
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

  %% linInterp
  rng(1);
  x = rand( 10, 1 );  x = sort( x );
  v = rand( 10, 1 );
  xq = rand(4,1);  xq = [ xq; x(2); max(x) ];
  y = rand( numel( xq ), 1 );

  Av = linInterp( x, v, xq );
  Av_correct = interp1( x, v, xq, 'linear', 0 );
  errInterp = norm( Av(:) - Av_correct(:) ) / norm( Av(:) );
  if errInterp > 1d-14
    error([ 'linInterp failed with error ', num2str( errInterp ) ]);
  end
  ATy = linInterp( x, y, xq, 'op', 'transp' );
  dp1 = dotP( Av, y );
  dp2 = dotP( v, ATy );
  err = abs( dp1 - dp2 ) / abs( dp1 );
  if err > 1d-14
    error([ 'linInterp adjoint failed with error ', num2str( err ) ]);
  else
    disp( 'linInterp passed' );
  end


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

  %% matchingPursuit
  close all
  fprintf( '\nTesting matchingPursuit: \n' );
  M = 30;  N = 20;
  K = 3;
  A = rand( M, N );
  nA = norms( A, 2, 1 );
  A = A ./ ( ones(M,1) * nA );
  x = zeros( N, 1 );
  x(1) = 1;
  x(8) = 3;
  x(N-2) = 2.10;
  b = A * x;
  xHat = matchingPursuit( A, b, K );
  err = norm( x(:) - xHat(:) ) / norm( x(:) );
  if err > 1
    error([ 'matchingPursuit failed with error ', num2str(err) ]);
  else
    disp( 'matchingPursuit passed' );
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

  %% orthogonalMatchingPursuit
  close all
  fprintf( '\nTesting orthogonalMatchingPursuit: \n' );
  M = 50;  N = 30;
  K = 3;
  A = rand( M, N );
  nA = norms( A, 2, 1 );
  A = A ./ ( ones(M,1) * nA );
  x = zeros( N, 1 );
  x(1) = 1;
  x(2) = 3;
  x(N-2) = 2.10;
  b = A * x;
  xHat = orthogonalMatchingPursuit( A, b, K );
  err = norm( x(:) - xHat(:) ) / norm( x(:) );
  if err > 1d-8
    error([ 'orthogonalMatchingPursuit failed with error ', num2str(err) ]);
  else
    disp( 'orthogonalMatchingPursuit passed' );
  end

  %% padData
  fprintf('\nTesting padData: \n');
  padded = padData( [1 2 3], 5 );
  err = norm( padded - [0 1 2 3 0] );
  if err>0, error(['padData failed with error ', num2str(err)]); end
  padded = padData( [1 2 3], 4 );
  err = norm( padded - [0 1 2 3] );
  if err>0, error(['padData failed with error ', num2str(err)]); end
  padded = padData( [1 2], 4 );
  err = norm( padded - [0 1 2 0] );
  if err>0, error(['padData failed with error ', num2str(err)]); end
  padded = padData( [1 2], 3 );
  err = norm( padded - [1 2 0] );
  if err>0, error(['padData failed with error ', num2str(err)]); end

  % test is cropData and padData are adjoints of each other
  x = rand(5,1);  y = rand(3,1);
  err = dotP( cropData(x,3), y ) - dotP( x, padData(y,5) );
  if err>0, error(['padData failed with error ', num2str(err)]); end
  x = rand(6,1); y = rand(4,1);
  err = dotP( cropData(x,4), y ) - dotP( x, padData(y,6) );
  if err>0, error(['padData failed with error ', num2str(err)]); end
  x = rand(5,1); y = rand(4,1);
  err = dotP( cropData(x,4), y ) - dotP( x, padData(y,5) );
  if err>0, error(['padData failed with error ', num2str(err)]); end
  x = rand(6,1); y = rand(5,1);
  err = dotP( cropData(x,5), y ) - dotP( x, padData(y,6) );
  if err>0, error(['padData failed with error ', num2str(err)]); end
  disp('padData passed');

  %% parforProgress
  N = 10;
  p = parforProgress( N );
  parfor n=1:N
    p.progress( n );
    pause(rand); % Replace with real code
  end
  p.clean;

  %% pdhg
  fprintf('\nTesting pdhg: \n');
  %minimize ||Ax - b||_2^2 + lambda ||x||
  A = rand( 20, 5 ) * 1;
  x = rand( 5, 1 ) * 1;
  x0 = rand( size( x ) );
  b = A * x;
  lambda = 0;
  f = @(x) lambda * norm( x, 1 );
  g = @(x) 0.5 * norm( x - b, 2 )^2;
  proxf = @(x,t) softThresh( x, t * lambda );
  proxgConj = @(x,gamma) proxConjL2Sq( x, gamma, 1, b );

  [~,s,~] = svd( A );
  normA = sqrt( s(1,1) );
  tau = 0.01 / s(1);
  [xStar,objValues,relDiffs] = pdhg( x0, proxf, proxgConj, tau, 'A', A, 'normA', normA, ...
    'N', 100000, 'f', f, 'g', g, 'printEvery', 1000, 'tol', 1d-99, 'verbose', true );   %#ok<ASGLU>
  mse = norm( x(:) - xStar(:) )^2 / numel( x );
  disp( [ x xStar x-xStar ] );
  fprintf([ '\npdhg ran to completion with MSE ', num2str(mse), '\n' ]);

  %% pdhgAdaptive
  fprintf('\nTesting pdhgAdaptive: \n');
  %minimize ||Ax - b||_2^2 + lambda ||x||
  A = rand( 7, 3 ) * 10;
  x = rand( 3, 1 ) * 1;
  x0 = rand( size( x ) );
  b = A * x;
  lambda = 10;
  f = @(x) lambda * norm( x, 1 );
  g = @(x) 0.5 * norm( x - b, 2 )^2;
  proxf = @(x,t) softThresh( x, t * lambda );
  proxgConj = @(x,gamma) proxConjL2Sq( x, gamma, 1, b );

  [~,s,~] = svd( A );
  normA = sqrt( s(1,1) );
  tau = 100;
  sigma = 100;
  [xStar,objValues] = pdhgAdaptive( x0, proxf, proxgConj, tau, 'sigma', sigma, 'A', A, ...
    'N', 1000, 'normA', normA, 'f', f, 'g', g, 'verbose', true );   %#ok<ASGLU>
  err = norm( x(:) - xStar(:) );
  fprintf('\npdhgAdaptive ran to completion\n');

  %% pdhgWLS
  fprintf('\nTesting pdhgWLS: \n');
  %minimize ||Ax - b||_2^2 + lambda ||x||
  A = rand( 7, 3 ) * 100;
  x = rand( 3, 1 ) * 1;
  x0 = rand( size( x ) );
  b = A * x;
  lambda = 1;
  f = @(x) lambda * norm( x, 1 );
  g = @(x) 0.5 * norm( x - b, 2 )^2;
  proxf = @(x,t) softThresh( x, t * lambda );
  proxgConj = @(x,gamma) proxConjL2Sq( x, gamma, 1, b );
  xStar = pdhgWLS( x0, proxf, proxgConj, 'A', A, 'N', 10000, 'f', f, 'g', g );   %#ok<NASGU>
  err = norm( x(:) - xStar(:) );
  fprintf('\npdhgWLS ran to completion\n');

  
  %% plotnice
  x = 1:10; y = 2*x;
  figure; subplot(4,1,1); plotnice(x);
  subplot(4,1,2); plotnice(x,y);
  subplot(4,1,3); plotnice(x, y, 'ro');
  subplot(4,1,4); plotnice( x, y, 'g+', 'MarkerSize', 20 );
  disp('Plotnice passed');

  %% polyFit2
  x = 1:10;  y = 1:10;
  [xs,ys] = meshgrid( x, y );
  xOrder = 2;  yOrder = 2;
  c = rand( xOrder+1, yOrder+1 );
  p = evaluatePoly2( c, xs, ys );
  cFit = polyFit2( xs, ys, p, xOrder, yOrder );
  err = norm( c(:) - cFit(:) );
  if err > 1d-6
    error([ 'polyFit2 failed with error ', num2str(err) ]);
  else
    disp('polyFit2 passed');
  end

  %% powerIteration
  fprintf('\nTesting powerIteration: \n');
  M = rand(3);
  est = powerIteration( M, false );
  normM = norm( M );
  err = abs( est - normM ) / normM;
  if err > 1d-5, error('Power Iteration failed'); end

  symmM = M + M';
  [est,flag] = powerIteration( symmM, true );  %#ok<ASGLU>
  normSymmM = norm( symmM );
  err = abs( est - normSymmM ) / normM;
  if err > 1d-5, error('powerIteration failed'); end

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
  if err > 1d-5, error('Power Iteration failed'); end
  disp('powerIteration passed');

  function out = applyMsymm( x )
    out = symmM * x;
  end
  x0 = rand( size(symmM,2), 1 );
  est = powerIteration( @applyMsymm, true, x0 );
  err = abs( est - normSymmM ) / normSymmM;
  if err > 1d-5, error('Power Iteration failed'); end

  disp('powerIteration passed');

  
  %% quadRoots
  a=2; b=0; c=-1;
  rs = quadRoots( a, b, c );
  rs = sort( rs );
  r2s = roots( [a b c ] );
  r2s = sort( r2s );
  err = norm( rs - r2s, 2 );
  if err < 1d-14
    disp('quadRoots passed');
  else
    error('quadRoots failed');
  end


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

  %% rotImg
  x = rand( 256 );
  y = rand( 256 );
  angle = 30 * pi/180;
  Ax = rotImg( x, angle );
  ATy = rotImg( y, angle, 'op', 'transp' );
  dp1 = dotP( Ax, y );
  dp2 = dotP( x, ATy );
  err = abs( dp1 - dp2 ) / abs( dp1 );
  if err < 1d-12
    disp('rotImg passed');
  else
    error('rotImg failed');
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

  %% smoothData
  x = rand(50);
  y = rand(50);
  Ax = smoothData( x, [7 5] );
  ATy = smoothData( y, [7 5], 'op', 'transp' );
  tmp1 = dotP( x, ATy );
  tmp2 = dotP( Ax, y );
  err = abs( tmp1 - tmp2 );
  if err < 1d-10
    disp( 'smoothData passed' );
  else
    error([ 'smoothData failed with an error of ', num2str(err) ]);
  end

  %% smoothImg
  x = rand(50);
  y = rand(50);
  Ax = smoothImg( x, [7 5] );
  ATy = smoothImg( y, [7 5], 'op', 'transp' );
  tmp1 = dotP( x, ATy );
  tmp2 = dotP( Ax, y );
  err = abs( tmp1 - tmp2 );
  if err < 1d-10
    disp( 'smoothImg passed' );
  else
    error([ 'smoothImg failed with an error of ', num2str(err) ]);
  end

  %% ufft2
  x = rand( 128, 128 );
  y = rand( 128, 128 );
  Fx = ufft2( x );
  FHy = uifft2( y );
  d1 = dotP( x, FHy );
  d2 = dotP( Fx, y );
  err = norm( d1(:) - d2(:) ) / norm( d1(:) );
  disp( err );
  if err < 1d-12
    disp( 'uifft2 is the adjoint of ufft2' );
  else
    error( [ 'ufft2 failed with an error of ', num2str( err ) ] );
  end

  %% vol2obj
  xs = 1 : 10;
  ys = 1 : 10;
  zs = 1 : 10;
  vol = bsxfun( @times, xs' * ys, reshape( zs, [1 1 10] ) );
  vol = vol / max( vol(:) );
  vol2obj( vol );

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


  %% withinPolygon
  polygon = [ 0 0 1; ...
              1 0 0];
  p = [ 0.5  0; ...
        0.5 -1];
  status = withinPolygon( polygon, p );
  err = max( abs( status - [ 1 0 ] ) );
  if err == 0
    disp( 'withinPolygon passed' );
  else
    error( 'withinPolygon failed' );
  end
  
  
  %% withinPolyhedron
  triangles = cat( 3, [0 0 0; 0 0 1; 1 0 0], ...
                      [0 0 0; 0 0 1; 0 1 0], ...
                      [0 0 0; 1 0 0; 0 1 0], ...
                      [0 0 1; 1 0 0; 0 1 0] );
  p = [  0.1  0.1 0.1; ...
           2    0   0; ...
         0.3  0.3 0.1; ...
           1    1   1; ...
           5    3  -1; ...
        0.05 0.05   0; ...
         0.2  0.1 0.3; ...
           1    0   0; ];
  status = withinPolyhedron( triangles, p );


end
