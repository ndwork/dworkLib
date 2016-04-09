
function [du,dv] = opticalFlow2D( img1, img2, varargin )
  % [du,dv] = opticalFlow2D( img1, img2 [, eta] )
  %
  % Computes the Horn Schunk optical flow vectors between images 1 and 2
  %
  % Inputs:
  % img1/img2 - 2D arrays
  %
  % Optional Inputs:
  % eta - smoothing parameter
  %
  % Outputs:
  % du - a 2D array representing the horiztonal components of the optical
  %   flow vectors
  % dv - a 2D array representing the vertical components of the optical
  %   flow vectors
  %
  % Written by Nicholas Dwork - Copyright 2016

  defaultEta = 1d-3;
  p = inputParser;
  p.addOptional( 'eta', defaultEta );
  p.parse( varargin{:} );
  eta = p.Results.eta;


  pyramid1 = makeImagePyramid( img1 );
  pyramid2 = makeImagePyramid( img2 );

  sLevel1 = size( pyramid1{end} );
  du = zeros( sLevel1(1), sLevel1(2) );
  dv = zeros( sLevel1(1), sLevel1(2) );

  for level=numel(pyramid1):-1:1
    disp(['Working on pyramid level ', num2str(level)]);
    tmp1 = pyramid1{level};
    [nRows, nCols] = size( pyramid1{level} );

    du = medfilt2( du, [3 3] );
    dv = medfilt2( dv, [3 3] );

    scale = nRows / size(du,1);
    du = imresize(du, [nRows nCols], 'bilinear') * scale;
    dv = imresize(dv, [nRows nCols], 'bilinear') * scale;

    tmp2 = ofInterp2D( pyramid2{level}, du, dv );

    [newDu, newDv] = ofADMM( tmp1, tmp2, eta );

    du = du + newDu;
    dv = dv + newDv;

    showDiagnostics = 0;
    if showDiagnostics==1
      close all;
      interped2 = ofInterp2D( pyramid2{level}, du, dv );
      figure, imshow( imresize( tmp1, ceil(512/nRows), 'nearest' ), [] );
      title('img1', 'FontSize', 20);
      figure, imshow( imresize( tmp2, ceil(512/nRows), 'nearest' ), [] );
      title('img2', 'FontSize', 20);
      figure, imshow( imresize( interped2, ceil(512/nRows), 'nearest' ), [] );
      title('interped2', 'FontSize', 20);
      figure, imshow( imresize( du, ceil(512/nRows), 'nearest' ), [] );
      title('du', 'FontSize', 20);
      figure, imshow( imresize( dv, ceil(512/nRows), 'nearest' ), [] );
      title('dv', 'FontSize', 20);
      drawnow;
    end

  end

end


function [du,dv] = ofADMM( img1, img2, eta )

  % ADMM Parameters
  rho = 0.001;

  [nRows, nCols] = size( img1 );

  applyD1 = @(u) cat(2, u(:,2:end) - u(:,1:end-1,:), zeros(nRows,1));
  applyD2 = @(u) cat(1, u(2:end,:) - u(1:end-1,:), zeros(1,nCols));
  applyD1Trans = @(u) cat(2, -u(:,1), u(:,1:end-2) - u(:,2:end-1), u(:,end-1));
  applyD2Trans = @(u) cat(1, -u(1,:), u(1:end-2,:) - u(2:end-1,:), u(end-1,:));

  applyM = @(u) applyD1Trans(applyD1(u)) + applyD2Trans(applyD2(u)) + u;
  allOnes = ones(nRows,nCols);
  eigValsM = dct2( applyM( idct2( allOnes ) ) );
  eigValsMInv = 1 ./ eigValsM;

  [Iu1, Iv1] = imgDeriv( img1 );
  [Iu2, Iv2] = imgDeriv( img2 );
  Iu = ( Iu1 + Iu2 ) / 2;
  Iv = ( Iv1 + Iv2 ) / 2;
  It = img2 - img1;

  % Initializations
  x1 = zeros( nRows, nCols );
  x2 = zeros( nRows, nCols );
  y1 = zeros( nRows, nCols );
  y2 = zeros( nRows, nCols );
  z1_1 = zeros( nRows, nCols );
  z1_2 = zeros( nRows, nCols );
  z2_1 = zeros( nRows, nCols );
  z2_2 = zeros( nRows, nCols );

  lambda1_1 = zeros( nRows, nCols );
  lambda1_2 = zeros( nRows, nCols );
  lambda2_1 = zeros( nRows, nCols );
  lambda2_2 = zeros( nRows, nCols );
  lambda3_1 = zeros( nRows, nCols );
  lambda3_2 = zeros( nRows, nCols );

  b = -It;
  Iub = Iu.*b;
  Ivb = Iv.*b;
  
  M11 = Iu.*Iu + rho;   M12 = Iu.*Iv;
  M21 = M12;            M22 = Iv.*Iv + rho;

  nIter = 500;
  objectives = zeros(nIter,1);
  for i=1:nIter

    objectives(i) = ofObjective( Iu, Iv, b, x1, x2, eta );

    % Update x
    arg1 = y1 + applyD1Trans(z1_1) + applyD2Trans(z1_2) - ...
      lambda1_1/rho - ...
      ( applyD1Trans(lambda2_1) + applyD2Trans(lambda2_2) )/rho;
    x1 = idct2( eigValsMInv .* dct2(arg1) );

    arg2 = y2 + applyD1Trans(z2_1) + applyD2Trans(z2_2) - ...
      lambda1_2/rho - ...
      ( applyD1Trans(lambda3_1) + applyD2Trans(lambda3_2) )/rho;
    x2 = idct2( eigValsMInv .* dct2(arg2) );

    % Update y
    nu1 = Iub + lambda1_1 + rho*x1;
    nu2 = Ivb + lambda1_2 + rho*x2;
    y2 = (nu2 - M21./M11.*nu1) ./ ( M22 - M21./M11.*M12 );
    y1 = 1./M11.*nu1 - 1./M11.*M12.*y2;

    % Update z
    z1_1 = softThresh( applyD1(x1) + lambda2_1/rho, eta/rho );
    z1_2 = softThresh( applyD2(x1) + lambda2_2/rho, eta/rho );
    z2_1 = softThresh( applyD1(x2) + lambda3_1/rho, eta/rho );
    z2_2 = softThresh( applyD2(x2) + lambda3_2/rho, eta/rho );

    % Update lambdas
    lambda1_1 = lambda1_1 + rho * ( x1 - y1 );
    lambda1_2 = lambda1_2 + rho * ( x2 - y2 );
    lambda2_1 = lambda2_1 + rho * ( applyD1(x1) - z1_1 );
    lambda2_2 = lambda2_2 + rho * ( applyD2(x1) - z1_2 );
    lambda3_1 = lambda3_1 + rho * ( applyD1(x2) - z2_1 );
    lambda3_2 = lambda3_2 + rho * ( applyD2(x2) - z2_2 );
  end

  du = reshape( x1, [nRows nCols] );
  dv = reshape( x2, [nRows nCols] );

  % for diagnostic purposes
  showDiagnostics = 0;
  if showDiagnostics==1
    close all;
    admmOptVal = ofObjective( Iu, Iv, b, x1, x2, eta );
    disp(['ADMM Optimal Value: ', num2str(admmOptVal) ] );
    ofRes = ofResidual( Iu, Iv, It, du, dv );
    figure, imshow( imresize( ofRes, ceil(512/nRows), 'nearest' ), [] );
    title('OF Residual', 'FontSize', 20);
    disp(['Min objective: ', num2str(min( objectives ))]);
    close all; figure; plot( objectives );
    %load( 'star.mat' );
    %objStar = ofObjective( Au, Av, b, x1Star, x2Star, eta );
    %figure, plot( ( objectives - objStar ) / objStar );
    %title('relative error vs iteration');
    drawnow;
  end

end


function out = augLagrange( Iu, Iv, It, x1, x2, y1, y2, ...
  z1_1, z1_2, z2_1, z2_2, eta, rho, ...
  lambda1_1, lambda1_2, lambda2_1, lambda2_2, lambda3_1, lambda3_2 )

  [nRows,nCols] = size( Iu );

  applyD1 = @(u) cat(2, u(:,2:end) - u(:,1:end-1,:), zeros(nRows,1));
  applyD2 = @(u) cat(1, u(2:end,:) - u(1:end-1,:), zeros(1,nCols));

  b = -It(:);
  Ay = Iu .* y1 + Iv .* y2;
  fOfY = 0.5*norm(Ay(:) - b, 2)^2;

  gOfZ = eta * norm(z1_1(:),1) + eta * norm(z1_2(:),1) + ...
    eta * norm(z2_1(:),1) + eta * norm(z2_2(:),1);

  tmpLam1_1 = lambda1_1 .* (x1-y1);
  tmpLam1_2 = lambda1_2 .* (x2-y2);
  costLam1 = sum( tmpLam1_1(:) ) + sum( tmpLam1_2(:) );

  tmpLam2 = lambda2_1 .* ( applyD1(x1) - z1_1 ) + ...
    lambda2_2 .* ( applyD2(x1) - z1_2 );
  costLam2 = sum( tmpLam2(:) );
  tmpLam3 = lambda3_1 .* ( applyD1(x2) - z2_1 ) + ...
    lambda3_2 .* ( applyD2(x2) - z2_2 );
  costLam3 = sum( tmpLam3(:) );

  D1x1 = applyD1(x1);   D1x2 = applyD1(x2);
  D2x1 = applyD2(x1);   D2x2 = applyD2(x2);
  augXY1 = norm( x1(:) - y1(:), 2 )^2;
  augXY2 = norm( x2(:) - y2(:), 2 )^2;
  augXZ1_1 = norm( D1x1(:) - z1_1(:), 2 )^2;
  augXZ1_2 = norm( D2x1(:) - z1_2(:), 2 )^2;
  augXZ2_1 = norm( D1x2(:) - z2_1(:), 2 )^2;
  augXZ2_2 = norm( D2x2(:) - z2_2(:), 2 )^2;
  
  out = fOfY + gOfZ + costLam1 + costLam2 + costLam3 + ...
    rho/2 * ( augXY1 + augXY2 + augXZ1_1 + augXZ1_2 + augXZ2_1 + augXZ2_2 );
  
  disp(['Aug Lag: ', num2str(out)]);
end

function [dx, dy] = imgDeriv( img )
  [M, N] = size( img );
  dx = zeros(M,N);
  dy = zeros(M,N);

  dx(:,2:N-1) = ( img(:,3:N) - img(:,1:N-2) ) / 2;
  dx(:,1) = img(:,2) - img(:,1);
  dx(:,N) = img(:,N) - img(:,N-1);

  dy(2:M-1,:) = ( img(3:M,:) - img(1:M-2,:) ) / 2;
  dy(1,:) = img(2,:) - img(1,:);
  dy(M,:) = img(M,:) - img(M-1,:);
end


function ofRes = ofResidual( Iu, Iv, It, du, dv )
  ofRes = Iu.*du + Iv.*dv + It;
end


function out = softThresh( in, thresh )
  out = sign(in) .* max( ( abs(in) - thresh ), 0 );
end


function out = ofObjective( Iu, Iv, b, du, dv, eta )
  [nRows,nCols] = size(Iu);

  applyD1 = @(u) cat(2, u(:,2:end) - u(:,1:end-1,:), zeros(nRows,1));
  applyD2 = @(u) cat(1, u(2:end,:) - u(1:end-1,:), zeros(1,nCols));

  D1du = applyD1(du);  D1dv = applyD1(dv);
  D2du = applyD2(du);  D2dv = applyD2(dv);
  
  Agam = Iu.*du + Iv.*dv;
  Agamb = Agam - b;
  out = 0.5*sum(Agamb(:).*Agamb(:)) + ...
    eta*norm(D1du(:),1) + eta*norm(D2du(:),1) + ...
    eta*norm(D1dv(:),1) + eta*norm(D2dv(:),1);
end


function pyramid = makeImagePyramid( img, nLevels, spacing )

  if nargin < 2
    nLevels = 4;
  end
  if nargin < 3
    spacing = 2;
  end

  smooth_sigma = spacing/2;
  f = fspecial('gaussian', 2*round(1.5*smooth_sigma) +1, smooth_sigma);
  ratio = 1. / spacing;

  pyramid = cell(nLevels,1);
  tmp = img;
  for m = 1:nLevels
    pyramid{m} = tmp;
    tmp = imfilter(tmp, f, 'corr', 'symmetric', 'same');  % Gauss filter
    tmp = imresize(tmp, ratio, 'bilinear');  % Downsampling
  end
end


