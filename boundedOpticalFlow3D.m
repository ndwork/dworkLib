
function [du,dv,dw] = boundedOpticalFlow3D( data1, data2, varargin )
  % [du,dv,dw] = boundedOpticalFlow3D( data1, data2, [, 'bound', bound, ...
  %   'eta', eta, 'rho', rho ] )
  %
  % Computes optical flow vectors
  %
  % Inputs:
  % data1/data2 - 3D arrays; optical flow vectors will point from voxels
  %   in data1 to voxels in data2
  %
  % Optional Inputs:
  % bound is a real positive number specifying the largest possible
  %   optical flow vector
  %
  % Written by Nicholas Dwork - Copyright 2015

  defaultBound = 1.2;
  defaultEta =  8d-10;
  defaultRho = 1d-3;
  p = inputParser;
  p.addRequired( 'data1' );
  p.addRequired( 'data2' );
  p.addParamValue( 'bound', defaultBound, @isnumeric );
  p.addParamValue( 'eta', defaultEta, @isnumeric );
  p.addParamValue( 'rho', defaultRho, @isnumeric );
  p.parse( data1, data2, varargin{:} );
  bound = p.Results.bound;
  eta = p.Results.eta;
  rho = p.Results.rho;

  pyramid1 = makeDataPyramid3( data1, 3 );
  pyramid2 = makeDataPyramid3( data2, 3 );

  sLevel1 = size( pyramid1{end} );
  du = zeros( sLevel1(1), sLevel1(2), sLevel1(3) );
  dv = zeros( sLevel1(1), sLevel1(2), sLevel1(3) );
  dw = zeros( sLevel1(1), sLevel1(2), sLevel1(3) );

  for level=numel(pyramid1):-1:1
    disp(['Working on pyramid level ', num2str(level)]);
    tmp1 = pyramid1{level};
    [nRows, nCols, nPages] = size( pyramid1{level} );

    du = medfilt3( du, [5 5 5], 'symmetric' );
    dv = medfilt3( dv, [5 5 5], 'symmetric' );
    dw = medfilt3( dw, [5 5 5], 'symmetric' );

    scale = nRows / size(du,1);
    du = resize(du, [nRows nCols nPages], 'bilinear') * scale;
    dv = resize(dv, [nRows nCols nPages], 'bilinear') * scale;
    dw = resize(dw, [nRows nCols nPages], 'bilinear') * scale;

    tmp2 = pyramid2{level};
    interp2 = ofInterp3D( tmp2, du, dv, dw );

    scaledEta = eta * numel(tmp1(:));
    [newDu,  newDv, newDw] = boundedOfADMM3D( tmp1, interp2, ...
      scaledEta, rho, bound );

    du = du + newDu;
    dv = dv + newDv;
    dw = dw + newDw;

    showDiagnostics = 1;
    if showDiagnostics==1
      close all;
      [nRows,nCols,nPages] = size(tmp1);
      midRow = ceil(nRows/2);
      midCol = ceil(nCols/2);
      midPage = ceil(nPages/2);
      img1 = squeeze(tmp1(:,:,midPage));
      minImg1 = min(img1(:));  maxImg1 = max(img1(:));
      figure;
      imshow( imresize( img1, ceil(512/nRows), 'nearest' ), ...
        [minImg1 maxImg1] );
      title('img1', 'FontSize', 20);
      img2 = squeeze(tmp2(:,:,midPage));
      figure;
      imshow( imresize( img2, ceil(512/nRows), 'nearest' ), ...
        [minImg1 maxImg1] );
      title('img2', 'FontSize', 20);
      interped2 = ofInterp3D( tmp2, du, dv, dw );
      interpedImg2 = squeeze(interped2(:,:,midPage));
      figure;
      imshow( imresize( interpedImg2, ceil(512/nRows), 'nearest' ), ...
        [minImg1 maxImg1] );
      title('interp2', 'FontSize', 20);
      duImg = squeeze(du(:,:,midPage));
      figure;
      imshow( imresize( duImg, ceil(512/nRows), 'nearest' ), [] );
      title('du', 'FontSize', 20);
      dvImg = squeeze(dv(:,:,midPage));
      figure;
      imshow( imresize( dvImg, ceil(512/nRows), 'nearest' ), [] );
      title('dv', 'FontSize', 20);
      dwImg = squeeze(dw(:,:,midPage));
      figure;
      imshow( imresize( dwImg, ceil(512/nRows), 'nearest' ), [] );
      title('dw', 'FontSize', 20);
      drawnow;
    end
    
  end

end


function [du,dv,dw] = boundedOfADMM3D( data1, data2, eta, rho, bound )

  % ADMM parameters
  nIter = 500;

  % data1 and data2 must be of type double
  data1 = double( data1 );
  data2 = double( data2 );

  % Check to make sure values of data1 and data2 are in a reasonable range
  if( max(data1(:)) > 1.5 || max(data2(:)) > 1.5 )
    error('Input values should be less than ~1');
  end

  [nRows, nCols, nPages] = size( data1 );

  applyD1 = @(u) cat(2, u(:,2:end,:) - u(:,1:end-1,:),zeros(nRows,1,nPages));
  applyD2 = @(u) cat(1, u(2:end,:,:) - u(1:end-1,:,:),zeros(1,nCols,nPages));
  applyD3 = @(u) cat(3, u(:,:,2:end) - u(:,:,1:end-1),zeros(nRows,nCols,1));
  applyD1Trans = @(u) cat(2,-u(:,1,:),u(:,1:end-2,:) - u(:,2:end-1,:),u(:,end-1,:));
  applyD2Trans = @(u) cat(1,-u(1,:,:),u(1:end-2,:,:) - u(2:end-1,:,:),u(end-1,:,:));
  applyD3Trans = @(u) cat(3,-u(:,:,1),u(:,:,1:end-2) - u(:,:,2:end-1),u(:,:,end-1));
  applyM = @(u) 2*u + applyD1Trans(applyD1(u)) + ...
    applyD2Trans(applyD2(u)) + applyD3Trans(applyD3(u));
  allOnes = ones(nRows,nCols,nPages);
  eigValsM = mirt_dctn(applyM(mirt_idctn(allOnes)));
  eigValsMInv = eigValsM.^(-1);

  [Iu1, Iv1, Iw1] = imgDeriv3D( data1 );
  [Iu2, Iv2, Iw2] = imgDeriv3D( data2 );
  Iu = 0.5 * ( Iu1 + Iu2 );
  Iv = 0.5 * ( Iv1 + Iv2 );
  Iw = 0.5 * ( Iw1 + Iw2 );
  It = data2 - data1;


  % Initializations
  x1 = zeros( nRows, nCols, nPages );
  x2 = zeros( nRows, nCols, nPages );
  x3 = zeros( nRows, nCols, nPages );
  r1 = zeros( nRows, nCols, nPages );
  r2 = zeros( nRows, nCols, nPages );
  r3 = zeros( nRows, nCols, nPages );
  y1 = zeros( nRows, nCols, nPages );
  y2 = zeros( nRows, nCols, nPages );
  y3 = zeros( nRows, nCols, nPages );
  z1_1 = zeros( nRows, nCols, nPages );
  z1_2 = zeros( nRows, nCols, nPages );
  z1_3 = zeros( nRows, nCols, nPages );
  z2_1 = zeros( nRows, nCols, nPages );
  z2_2 = zeros( nRows, nCols, nPages );
  z2_3 = zeros( nRows, nCols, nPages );
  z3_1 = zeros( nRows, nCols, nPages );
  z3_2 = zeros( nRows, nCols, nPages );
  z3_3 = zeros( nRows, nCols, nPages );

  lambda1u = zeros( nRows, nCols, nPages );
  lambda1v = zeros( nRows, nCols, nPages );
  lambda1w = zeros( nRows, nCols, nPages );
  lambda2u = zeros( nRows, nCols, nPages );
  lambda2v = zeros( nRows, nCols, nPages );
  lambda2w = zeros( nRows, nCols, nPages );
  lambda3u_1 = zeros( nRows, nCols, nPages );
  lambda3u_2 = zeros( nRows, nCols, nPages );
  lambda3u_3 = zeros( nRows, nCols, nPages );
  lambda3v_1 = zeros( nRows, nCols, nPages );
  lambda3v_2 = zeros( nRows, nCols, nPages );
  lambda3v_3 = zeros( nRows, nCols, nPages );
  lambda3w_1 = zeros( nRows, nCols, nPages );
  lambda3w_2 = zeros( nRows, nCols, nPages );
  lambda3w_3 = zeros( nRows, nCols, nPages );

  b = -It;
  Iub = Iu.*b;
  Ivb = Iv.*b;
  Iwb = Iw.*b;

  M11 = Iu.*Iu + rho;   M12 = Iu.*Iv;           M13 = Iu.*Iw;
  M21 = M12;            M22 = Iv.*Iv + rho;     M23 = Iv.*Iw;
  M31 = M13;            M32 = M23;              M33 = Iw.*Iw + rho;

  KInv = ( M22 - M21./M11.*M12 );
  K = 1. / KInv;
  CInv = -M31./M11.*M12.*K.*M21./M11.*M13 + M31./M11.*M12.*K.*M23 - ...
    M31./M11.*M13 + M32.*K.*M21./M11.*M13 - M32.*K.*M23 + M33;
  C = 1. / CInv;

  objectives = zeros(nIter,1);
  for i=1:nIter
    if mod(i,100)==0
      disp([ 'Working on iteration ', num2str(i), ' of ', num2str(nIter) ]);
    end

    thisObjective = ofObjective( Iu, Iv, Iw, b, x1, x2, x3, eta );
    objectives(i) = thisObjective;
    
    % Update x
    arg1 = r1 + y1 + ...
      applyD1Trans(z1_1) + applyD2Trans(z1_2) + applyD3Trans(z1_3) ...
      - ( lambda1u + lambda2u ) / rho ...
      - ( applyD1Trans(lambda3u_1) + applyD2Trans(lambda3u_2) ...
        + applyD3Trans(lambda3u_3) ) / rho;
    x1 = mirt_idctn(eigValsMInv.*mirt_dctn(arg1));
    
    arg2 = r2 + y2 + ...
      applyD1Trans(z2_1) + applyD2Trans(z2_2) + applyD3Trans(z2_3) ...
      - ( lambda1v + lambda2v ) / rho ...
      - ( applyD1Trans(lambda3v_1) + applyD2Trans(lambda3v_2) ...
        + applyD3Trans(lambda3v_3) ) / rho;
    x2 = mirt_idctn(eigValsMInv.*mirt_dctn(arg2));
    
    arg3 = r3 + y3 + ...
      applyD1Trans(z3_1) + applyD2Trans(z3_2) + applyD3Trans(z3_3) ...
      - ( lambda1w + lambda2w ) / rho ...
      - ( applyD1Trans(lambda3w_1) + applyD2Trans(lambda3w_2) ...
        + applyD3Trans(lambda3w_3) ) / rho;
    x3 = mirt_idctn(eigValsMInv.*mirt_dctn(arg3));

    % Update r
    rArg1 = x1 + lambda1u / rho;
    rArg2 = x2 + lambda1v / rho;
    rArg3 = x3 + lambda1w / rho;
    rMag = sqrt( rArg1.*rArg1 + rArg2.*rArg2 + rArg3.*rArg3 );
    factor = min( bound ./ rMag, 1 );
    r1 = rArg1 .* factor;
    r2 = rArg2 .* factor;
    r3 = rArg3 .* factor;
    
    % Update y
    nu1 = Iub + lambda2u + rho*x1;
    nu2 = Ivb + lambda2v + rho*x2;
    nu3 = Iwb + lambda2w + rho*x3;
    y3 = C.*( nu3 - M31./M11.*nu1 + M31./M11.*M12.*K.*nu2 - ...
      M31./M11.*M12.*K.*M21./M11.*nu1 - M32.*K.*nu2 + ...
      M32.*K.*M21./M11.*nu1 );
    y2 = K .* (nu2 - M21./M11.*nu1 + M21./M11.*M13.*y3 - M23.*y3 );
    y1 = 1./M11 .* ( nu1 - M12.*y2 - M13.*y3 );
    
    % Update z
    z1_1 = softThresh( applyD1(x1) + lambda3u_1/rho, eta/rho );
    z1_2 = softThresh( applyD2(x1) + lambda3u_2/rho, eta/rho );
    z1_3 = softThresh( applyD3(x1) + lambda3u_3/rho, eta/rho );
    z2_1 = softThresh( applyD1(x2) + lambda3v_1/rho, eta/rho );
    z2_2 = softThresh( applyD2(x2) + lambda3v_2/rho, eta/rho );
    z2_3 = softThresh( applyD3(x2) + lambda3v_3/rho, eta/rho );
    z3_1 = softThresh( applyD1(x3) + lambda3w_1/rho, eta/rho );
    z3_2 = softThresh( applyD2(x3) + lambda3w_2/rho, eta/rho );
    z3_3 = softThresh( applyD3(x3) + lambda3w_3/rho, eta/rho );
    
    % Update lambdas
    lambda1u = lambda1u + rho * ( x1 - r1 );
    lambda1v = lambda1u + rho * ( x2 - r2 );
    lambda1w = lambda1u + rho * ( x3 - r3 );
    lambda2u = lambda2u + rho * ( x1 - y1 );
    lambda2v = lambda2v + rho * ( x2 - y2 );
    lambda2w = lambda2w + rho * ( x3 - y3 );
    lambda3u_1 = lambda3u_1 + rho * ( applyD1(x1) - z1_1 );
    lambda3u_2 = lambda3u_2 + rho * ( applyD2(x1) - z1_2 );
    lambda3u_3 = lambda3u_3 + rho * ( applyD3(x1) - z1_3 );
    lambda3v_1 = lambda3v_1 + rho * ( applyD1(x2) - z2_1 );
    lambda3v_2 = lambda3v_2 + rho * ( applyD2(x2) - z2_2 );
    lambda3v_3 = lambda3v_3 + rho * ( applyD3(x2) - z2_3 );
    lambda3w_1 = lambda3w_1 + rho * ( applyD1(x3) - z3_1 );
    lambda3w_2 = lambda3w_2 + rho * ( applyD2(x3) - z3_2 );
    lambda3w_3 = lambda3w_3 + rho * ( applyD3(x3) - z3_3 );

  end

  du = x1;  dv = x2;  dw = x3;

  % for diagnostics
  showDiagnostics = 0;
  if showDiagnostics==1
    close all;
    admmOptVal = ofObjective( Iu, Iv, Iw, b, x1, x2, x3, eta );
    disp(['ADMM Optimal Value: ', num2str(admmOptVal) ] );
    figure; plot( objectives );
    ofRes = ofResidual3D( Iu, Iv, Iw, It, du, dv, dw );
    %figure, imshow( imresize( ofRes, ceil(512/nRows), 'nearest' ), [] );
    %title('OF Residual', 'FontSize', 20);
    %load( 'star3D.mat' );
    %objStar = ofObjective( Iu, Iv, Iw, b, x1Star, x2Star, x3Star, ...
    %  eta, D1, D2, D3 );
    %figure, plot( ( objectives - objStar ) / objStar );
    %title('relative error vs iteration');
    drawnow;

    disp(['Rho: ', num2str(rho)]);
    disp(['Max objective: ', num2str(max( objectives ))]);
    disp(['Min objective: ', num2str(min( objectives ))]);
    close all; figure; plot( objectives ); drawnow;
  end

end


function out = augLagrange( Iu, Iv, Iw, It, x1, x2, x3, r1, r2, r3, ...
  y1, y2, y3, z1_1, z1_2, z1_3, z2_1, z2_2, z2_3, z3_1, z3_2, z3_3, ...
  eta, rho, bound, ...
  lambda1u, lambda1v, lambda1w, lambda2u, lambda2v, lambda2w, ...
  lambda3u_1, lambda3u_2, lambda3u_3, lambda3v_1, lambda3v_2, ...
  lambda3v_3, lambda3w_1, lambda3w_2, lambda3w_3 )

  rMag = sqrt( r1.*r1 + r2.*r2 + r3.*r3 );
  maxRMag = max( rMag(:) );
  if maxRMag > bound
    out = Inf;
    return
  end

  [nRows,nCols,nPages] = size( Iu );

  applyD1 = @(u) cat(2, u(:,2:end,:) - u(:,1:end-1,:), zeros(nRows,1,nPages));
  applyD2 = @(u) cat(1, u(2:end,:,:) - u(1:end-1,:,:), zeros(1,nCols,nPages));
  applyD3 = @(u) cat(3, u(:,:,2:end) - u(:,:,1:end-1), zeros(nRows,nCols,1));

  D1x1 = applyD1(x1);   D1x2 = applyD1(x2);   D1x3 = applyD1(x3);
  D2x1 = applyD2(x1);   D2x2 = applyD2(x2);   D2x3 = applyD2(x3);
  D3x1 = applyD3(x1);   D3x2 = applyD3(x2);   D3x3 = applyD3(x3);
  
  b = -It(:);
  Ay = Iu .* y1 + Iv .* y2 + Iw .* y3;
  gOfW = 0.5*norm(Ay(:) - b, 2)^2 + ...
    eta*norm(z1_1(:),1) + eta*norm(z1_2(:),1) + eta*norm(z1_3(:),1) + ...
    eta*norm(z2_1(:),1) + eta*norm(z2_2(:),1) + eta*norm(z2_3(:),1) + ...
    eta*norm(z3_1(:),1) + eta*norm(z3_2(:),1) + eta*norm(z3_3(:),1);

  tmpLam1u = lambda1u .* (x1 - r1);
  tmpLam1v = lambda1v .* (x2 - r2);
  tmpLam1w = lambda1w .* (x3 - r3);
  costLam1 = sum( tmpLam1u(:) + tmpLam1v(:) + tmpLam1w(:) );

  tmpLam2u = lambda2u .* (x1 - y1);
  tmpLam2v = lambda2v .* (x2 - y2);
  tmpLam2w = lambda2w .* (x3 - y3);
  costLam2 = sum( tmpLam2u(:) + tmpLam2v(:) + tmpLam2w(:) );

  tmpLam3u = lambda3u_1 .* ( D1x1 - z1_1 ) + ...
    lambda3u_2 .* ( D2x1 - z1_2 ) + ...
    lambda3u_3 .* ( D3x1 - z1_3 );
  costLam3u = sum( tmpLam3u(:) );

  tmpLam3v = lambda3v_1 .* ( D1x2 - z2_1 ) + ...
    lambda3v_2 .* ( D2x2 - z2_2 ) + ...
    lambda3v_3 .* ( D3x2 - z2_3 );
  costLam3v = sum( tmpLam3v(:) );

  tmpLam3w = lambda3w_1 .* ( D1x3 - z3_1 ) + ...
    lambda3w_2 .* ( D2x3 - z3_2 ) + ...
    lambda3w_3 .* ( D3x3 - z3_3 );
  costLam3w = sum( tmpLam3w(:) );

  augXR1 = norm( x1(:) - r1(:), 2 )^2;
  augXR2 = norm( x2(:) - r2(:), 2 )^2;
  augXR3 = norm( x3(:) - r3(:), 2 )^2;
  augXY1 = norm( x1(:) - y1(:), 2 )^2;
  augXY2 = norm( x2(:) - y2(:), 2 )^2;
  augXY3 = norm( x3(:) - y3(:), 2 )^2;
  augXZ1_1 = norm( D1x1(:) - z1_1(:), 2 )^2;
  augXZ1_2 = norm( D2x1(:) - z1_2(:), 2 )^2;
  augXZ1_3 = norm( D3x1(:) - z1_3(:), 2 )^2;
  augXZ2_1 = norm( D1x2(:) - z2_1(:), 2 )^2;
  augXZ2_2 = norm( D2x2(:) - z2_2(:), 2 )^2;
  augXZ2_3 = norm( D3x2(:) - z2_3(:), 2 )^2;
  augXZ3_1 = norm( D1x3(:) - z3_1(:), 2 )^2;
  augXZ3_2 = norm( D2x3(:) - z3_2(:), 2 )^2;
  augXZ3_3 = norm( D3x3(:) - z3_3(:), 2 )^2;

  out = gOfW + costLam1 + costLam2 + ...
    costLam3u + costLam3v + costLam3w + ...
    rho/2 * ( augXR1 + augXR2 + augXR3 + augXY1 + augXY2 + augXY3 + ...
      augXZ1_1 + augXZ1_2 + augXZ1_3 + augXZ2_1 + augXZ2_2 + augXZ2_3 + ...
      augXZ3_1 + augXZ3_2 + augXZ3_3 );

  disp(['Aug Lag: ', num2str(out)]);
end


function [dx, dy, dz] = imgDeriv3D( data )
  [M, N, K] = size( data );
  dx = zeros(M,N,K);
  dy = zeros(M,N,K);
  dz = zeros(M,N,K);

  dx(:,2:N-1,:) = ( data(:,3:N,:) - data(:,1:N-2,:) ) / 2;
  dx(:,1,:) = data(:,2,:) - data(:,1,:);
  dx(:,N,:) = data(:,N,:) - data(:,N-1,:);
  
  dy(2:M-1,:,:) = ( data(3:M,:,:) - data(1:M-2,:,:) ) / 2;
  dy(1,:,:) = data(2,:,:) - data(1,:,:);
  dy(M,:,:) = data(M,:,:) - data(M-1,:,:);
  
  dz(:,:,2:K-1) = ( data(:,:,3:K) - data(:,:,1:K-2) ) / 2;
  dz(:,:,1) = data(:,:,2) - data(:,:,1);
  dz(:,:,K) = data(:,:,K) - data(:,:,K-1);
end


function ofRes = ofResidual3D( Iu, Iv, Iw, It, du, dv, dw )
  ofRes = Iu.*du + Iv.*dv + It + Iw.*dw;
end


function out = softThresh( in, thresh )
  out = sign(in) .* max( ( abs(in) - thresh ), 0 );
end


function out = ofObjective( Iu, Iv, Iw, b, du, dv, dw, eta )
  [nRows,nCols,nPages] = size(Iu);

  applyD1 = @(u) cat(2, u(:,2:end,:) - u(:,1:end-1,:),zeros(nRows,1,nPages));
  applyD2 = @(u) cat(1, u(2:end,:,:) - u(1:end-1,:,:),zeros(1,nCols,nPages));
  applyD3 = @(u) cat(3, u(:,:,2:end) - u(:,:,1:end-1),zeros(nRows,nCols,1));

  D1du = applyD1(du);   D1dv = applyD1(dv);   D1dw = applyD1(dw);
  D2du = applyD2(du);   D2dv = applyD2(dv);   D2dw = applyD2(dw);
  D3du = applyD3(du);   D3dv = applyD3(dv);   D3dw = applyD3(dw);

  Agam = Iu.*du + Iv.*dv + Iw.*dw;
  Agamb = Agam - b;
  out = 0.5*sum(Agamb(:).*Agamb(:)) + ...
    eta*norm(D1du(:),1) + eta*norm(D2du(:),1) + eta*norm(D3du(:),1) + ...
    eta*norm(D1dv(:),1) + eta*norm(D2dv(:),1) + eta*norm(D3dv(:),1) + ...
    eta*norm(D1dw(:),1) + eta*norm(D2dw(:),1) + eta*norm(D3dw(:),1);
end


function pyramid = makeDataPyramid3( data, nLevels, spacing )
  if nargin < 2
    nLevels = 4;
  end
  if nargin < 3
    spacing = 2;
  end

  smooth_sigma = spacing/2;
  f = gaussFilt3( smooth_sigma, 2*round(1.5*smooth_sigma)+1 );
  ratio = 1. / spacing;

  pyramid = cell(nLevels,1);
  tmp = data;
  for m = 1:nLevels
    pyramid{m} = tmp;
    sTmp = size( tmp );
    newSize = floor( sTmp .* ratio );
    tmp = imfilter(tmp, f, 'corr', 'symmetric', 'same');  % Gauss filter
    tmp = resize(tmp, newSize, 'bilinear');  % Downsampling
  end
end


function h = gaussFilt3( sig, hsize )
  siz   = (hsize-1)/2;
  if ndims(siz) ~= 3 siz = ones(1,3)*siz; end;
  if ndims(sig) ~= 3 sig = ones(1,3)*sig; end;
  [x,y,z] = ndgrid(-siz(1):siz(1),-siz(2):siz(2),-siz(3):siz(3));
  h = exp(-(x.*x/2/sig(1)^2 + y.*y/2/sig(2)^2 + z.*z/2/sig(3)^2));
  h = h/sum(h(:));
end

