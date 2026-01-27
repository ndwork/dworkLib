
function [ K, Hs ] = calibrateCamFromSquares( squares, varargin )
  % K = calibrateCamFromSquares( squares )
  %
  % Written according to (8.12) from "Multiple View Geometry" by Hartley and Zisserman,
  % 2nd edition (page 211).
  % Note that it is suggested that the points of squares be normalized using normalizePts2D.
  %
  % For the case where zero skew is enforced, the method comes from Chapter 2 of "Emerging
  % Methods in Computer Vision" by Gérard Medioni and Sing Bing Kang.
  %
  % Inputs:
  % squares - a 3D array of size 4 x 2 x nSquares (nSquares must be >=3 ) that represents
  %   the image points of each square's corners.  The first/second column represents the
  %   first/second coordinate, respectively.  The points are projected onto [ (0,0),
  %   (1,0), (0,1), (1,1) ], so they should be in a comparable order.
  %
  % Optional Inputs:
  % zeroSkew - if true, then zero skew is assumed
  % principalPoint - a two element array.  If provided, then this principal point is assumed.
  %
  % Outputs:
  % K - the 3 x 3 intrisic camera matrix
  %
  % Written by Nicholas Dwork, Copyright 2023
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addParameter( 'doCheck', false );
  p.addParameter( 'zeroSkew', true );
  p.addParameter( 'principalPoint', [], @isnumeric );
  p.addParameter( 'tol', 1d-8, @isnonnegative );
  p.parse( varargin{:} );
  doCheck = p.Results.doCheck;
  zeroSkew = p.Results.zeroSkew;
  principalPoint = p.Results.principalPoint;
  tol = p.Results.tol;

  nSquares = size( squares, 3 );
  if size(squares,1) ~= 4 || size(squares,2) ~= 2 || nSquares < 3
    error('Input "squares" must be a 4 x 2 x nSquares array with nSquares >= 3.');
  end

  simplePts = [ 0 0; 1 0; 0 1; 1 1; ];
  Hs = zeros( 3, 3, nSquares );

  for squareIndx = 1 : nSquares
    Hs(:,:,squareIndx) = homographyFromPts2D( simplePts, squares(:,:,squareIndx) );

    if doCheck == true
      errH = checkH( Hs(:,:,squareIndx), squares(:,:,squareIndx) );
      if errH > 1d-11
        error([ 'Homography check for square ', num2str(squareIndx), ' failed.' ]);
      else
        disp([ 'Homography check for square ', num2str(squareIndx), ' passed.' ]);
      end
    end
  end

  if zeroSkew
    if numel( principalPoint ) > 0
      iac = findIAC_zeroSkew_withPP( Hs, principalPoint(1), principalPoint(2) );
    else
      iac = findIAC_zeroSkew( Hs );
    end
  else
    if numel( principalPoint ) > 0
      iac = findIAC_nonzeroSkew_wPP( Hs, principalPoint(1), principalPoint(2) );
    else
      iac = findIAC_nonzeroSkew( Hs );
    end
  end

  [~, D] = eig( iac );
  if any(diag(D) < 0)
    iac = -iac;
  end

  invIAC = inv( iac );

  if min([ invIAC(1,1), invIAC(2,2), invIAC(3,3) ]) < 0
    error( 'Inverse IAC is not positive definite.  There are numerical issues' );
  end

  try
    [ K, cholFlag ] = chol( invIAC, 'upper' );   %#ok<ASGLU>
  catch
    error('Inverse IAC is not positive definite.');
  end

  K = K ./ K(3,3);

  if zeroSkew
    if abs( K(1,2) ) > tol
      warning( 'large skew suggests numerical issues' );
    end
    K(1,2) = 0;
  end
end


%-- Support routines

function errH = checkH( H, squarePts )
  simplePts = [ 0 0; 1 0; 0 1; 1 1; ];
  Hsimple = H * [ simplePts'; 1 1 1 1; ];
  Hpts = Hsimple ./ Hsimple(3,:);
  Hpts = Hpts(1:2,:);
  diffPts = squarePts - Hpts';
  errH = norm( diffPts(:), 2 );
end


function iac = findIAC_nonzeroSkew( Hs )
  nSquares = size( Hs, 3 );
  A = zeros( 2*nSquares, 6 );

  for i = 1 : nSquares
    H = Hs(:,:,i);

    A(i,:) = [ H(1,1) * H(1,2) ...
               H(1,1) * H(2,2) + H(2,1) * H(1,2) ...
               H(1,1) * H(3,2) + H(3,1) * H(1,2) ...
               H(2,1) * H(2,2) ...
               H(2,1) * H(3,2) + H(3,1) * H(2,2) ...
               H(3,1) * H(3,2) ];

    d = H(:,1) - H(:,2);  % difference
    s = H(:,1) + H(:,2);  % sum
    A(i+nSquares,:) = [ d(1) * s(1)  ...
                        d(1) * s(2) + d(2) * s(1) ...
                        d(1) * s(3) + d(3) * s(1) ...
                        d(2) * s(2) ...
                        d(2) * s(3) + d(3) * s(2) ...
                        d(3) * s(3) ];
  end

  [~,s,v] = svd( A, 'vector' );   %#ok<ASGLU>

  iacVec = v(:,end);
  iac = [ iacVec(1) iacVec(2) iacVec(3); ...
          iacVec(2) iacVec(4) iacVec(5); ...
          iacVec(3) iacVec(5) iacVec(6); ];
end


function iac = findIAC_nonzeroSkew_wPP( Hs, px, py )
  n = size(Hs, 3);
  A = zeros(2*n + 2, 6);   % 2 constraints per square + 2 for principal point
  b = zeros(2*n + 2, 1);

  for i = 1:n
      H = Hs(:,:,i);

      % --- First constraint: h1^T ω h2 = 0 ---
      A(i, :) = [ H(1,1)*H(1,2), ...
                  H(1,1)*H(2,2) + H(2,1)*H(1,2), ...
                  H(1,1)*H(3,2) + H(3,1)*H(1,2), ...
                  H(2,1)*H(2,2), ...
                  H(2,1)*H(3,2) + H(3,1)*H(2,2), ...
                  H(3,1)*H(3,2) ];

      % --- Second constraint: h1^T ω h1 = h2^T ω h2 ---
      d = H(:,1) - H(:,2);   % difference
      s = H(:,1) + H(:,2);   % sum
      A(i+n, :) = [ d(1)*s(1), ...
                    d(1)*s(2) + d(2)*s(1), ...
                    d(1)*s(3) + d(3)*s(1), ...
                    d(2)*s(2), ...
                    d(2)*s(3) + d(3)*s(2), ...
                    d(3)*s(3) ];
  end

  % --- Principal point constraints ---
  % ω13 + px * ω11 = 0  →  px*ω11 + 1*ω13 = 0
  % ω23 + py * ω22 = 0  →  py*ω22 + 1*ω23 = 0

  A(2*n+1, :) = [ px,  0,  1,  0,  0,  0 ];   % px*ω11 + ω13 = 0
  A(2*n+2, :) = [  0,  0,  0, py,  1,  0 ];   % py*ω22 + ω23 = 0

  v = A \ b;   % 6×1 solution

  iac = [ v(1)  v(2)  v(3); ...
          v(2)  v(4)  v(5); ...
          v(3)  v(5)  v(6) ];
end


function iac = findIAC_zeroSkew( Hs )
  nSquares = size( Hs, 3 );
  A = zeros( 2*nSquares, 5 );
  for i = 1 : nSquares
      H = Hs(:,:,i);
      A(i,:) = [ H(1,1) * H(1,2) ...
                 H(1,1) * H(3,2) + H(3,1) * H(1,2) ...
                 H(2,1) * H(2,2) ...
                 H(2,1) * H(3,2) + H(3,1) * H(2,2) ...
                 H(3,1) * H(3,2) ];
      d = H(:,1) - H(:,2);
      s = H(:,1) + H(:,2);
      A(i+nSquares,:) = [ d(1) * s(1) ...
                          d(1) * s(3) + d(3) * s(1) ...
                          d(2) * s(2) ...
                          d(2) * s(3) + d(3) * s(2) ...
                          d(3) * s(3) ];
  end

  [~,s,v] = svd( A, 'vector' );   %#ok<ASGLU>

  iacVec = v(:,end);
  iac = [ iacVec(1)     0       iacVec(2); ...
             0       iacVec(3)  iacVec(4); ...
          iacVec(2)  iacVec(4)  iacVec(5); ];
end


function iac = findIAC_zeroSkew_withPP(Hs, px, py)
    n = size(Hs,3);
    A = zeros(2*n + 2, 5);
    b = zeros(2*n + 2, 1);

    for i=1:n
        H = Hs(:,:,i);
        A(i,:) = [ H(1,1)*H(1,2), ...
                   H(1,1)*H(3,2)+H(3,1)*H(1,2), ...
                   H(2,1)*H(2,2), ...
                   H(2,1)*H(3,2)+H(3,1)*H(2,2), ...
                   H(3,1)*H(3,2) ];
        d = H(:,1)-H(:,2); s = H(:,1)+H(:,2);
        A(i+n,:) = [ d(1)*s(1), d(1)*s(3)+d(3)*s(1), ...
                     d(2)*s(2), d(2)*s(3)+d(3)*s(2), d(3)*s(3) ];
    end

    % ---- principal-point constraints ----
    % ω13 + px*ω11 = 0   →  row:  px  1  0  0  0
    % ω23 + py*ω22 = 0   →  row:   0  0 py  1  0
    A(2*n+1,:) = [ px  1  0  0  0 ];
    A(2*n+2,:) = [  0  0 py  1  0 ];

    v = A \ b;
    iac = [ v(1)  0   v(2); ...
            0    v(3) v(4); ...
            v(2) v(4) v(5); ];
end
