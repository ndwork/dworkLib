function out = mri_reconSPIRiT( kData, kernel_sz, acr_sz, varargin )
  %
  % *********************
  %   Input Parameters:
  % *********************
  %
  %    kData:  A 3D (size ny x nx x ncoils) array of measured k-space. The
  %    assumption in GRAPPA is that there is a fully sampled region in the
  %    center of the image. The center point convention is as follows:
  %
  %    if our array is odd size, we choose the center point
  %           o o o x o o o
  %    if our array is even size, we choose length / 2 + 1
  %           o o o x o o 
  %
  %    kernel_sz: a 2 element vector containing the dimensions of the kernel.
  %    We use MATLAB dimensions, so it's [row column]
  %    For example, a 3x1 kernel will fit inside a 5x5 ACR as:
  %                 x o o o o
  %                 x o o o o
  %                 x o o o o
  %                 o o o o o 
  %                 o o o o o 
  %    acr_sz: A 2 element vector containing the dimensions of the ACR. ACR
  %    must be odd sized in both directions.
  %
  % Written by Alex McManus - Copyright 2022
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It is offered without any
  % warranty expressed or implied, including the implied warranties of merchantability or fitness
  % for a particular purpose.

  p = inputParser;
  p.addParameter( 'alg', 'fista_wLS', @(x) true );
  p.addParameter( 'checkProx', false );
  p.addParameter( 'eps', 0, @(x) numel(x) == 0 || isnonnegative(x) );
  p.addParameter( 'verbose', false );
  p.addParameter( 'weights', [], @isnumeric );
  p.addParameter( 'x0', [], @isnumeric );
  p.parse( varargin{:} );
  alg = p.Results.alg;
  eps = p.Results.eps;
  checkProx = p.Results.checkProx;
  verbose = p.Results.verbose;
  weights = p.Results.weights;
  x0 = p.Results.x0;

  if numel( kernel_sz ) == 1, kernel_sz = [ kernel_sz kernel_sz ]; end
  if numel( eps ) == 0, eps = 0; end

  nCoils = size( kData, 3 );
  
  if numel( weights ) == 0
    acr = cropData( kData, [ acr_sz nCoils ] );
    weights = mri_reconSPIRiT_get_weights( acr, kernel_sz );
  end

  idx_acq = kData~=0;

  y_data = kData(idx_acq);
  %eps = 1e-2;

  %x0 = mri_reconGRAPPA( kData, kernel_sz, acr_sz );
  if numel( x0 ) == 0, x0 = kData; end

  if checkProx == true
    proxTest = @(in, sc) projectOntoBall( in, sqrt(eps) );

    test1 = proxh(x0, 0);
    test2 = proxAffine( proxTest, x0, @applyD, -1*y_data, 1 );
    
    err = norm(test1(:) - test2(:))/norm(test2(:));
    if abs(err) > 1e-6
      error('Something wrong with proximal operator');
    end
  end

  if strcmp( alg, 'fista' )
    [xStar, objVal, relDiffs] = fista(0*x0(:), @grad_g, @proxh, 'g', @normG, 'h', @h, 'N', 100, ...
      'verbose', verbose, 't', 0.05 );   %#ok<ASGLU>

  elseif strcmp( alg, 'fista_wLS' )
    [xStar, objVal, relDiffs] = fista_wLS( x0(:), @g, @grad_g, @proxh, 'h', @h, 'N', 100, ...
      'verbose', verbose, 't0', 0.01);   %#ok<ASGLU>
  end

  out = reshape( xStar, size(kData) );

  function out = proxh( in, sc )   %#ok<INUSD>
    din = in( idx_acq ) - y_data;
    tmp_proj = projectOntoBall( din, sqrt(eps) );
    tmp2 = din - tmp_proj;
    out = zeros( size( kData ) );
    out( idx_acq ) = tmp2;
    out = in(:) - out(:);
  end

  function out = applyG( in, op )
    in = reshape( in, size( kData ) );
    if nargin < 2 || strcmp(op, 'notransp')
      out = spirit_conv( in, weights );
    else
      out = spirit_conv_adj( in, weights );
    end
    out = out(:);
  end

  function out = h(x)
    out = normIndicator( x, y_data, sqrt(eps) );
  end

  function out = normIndicator( x, y, r )
    dx = x( idx_acq );
    if norm(dx - y) < r
      out = 0;
    else
      out = Inf;
    end
  end

  % function out = projB( x, y, r )
  %   dxy = x - y;
  %   if norm(dxy) < r
  %     out = x;
  %   else
  %     out = dxy * (r / norm(dxy)) + y;
  %   end
  % end

  function out = g( in )
    tmpg = applyG( in );
    out = 0.5 * norm( tmpg(:) - in(:) )^2;
  end

  function out = grad_g(in)
    tmpgrad = applyG( in ) - in;
    out = applyG( tmpgrad, 'transp' ) - tmpgrad;
    out = out(:);
  end

  function out = applyD( in, op )  % Apply the sampling mask
    if strcmp( op, 'transp' )
      out = 0 * kData;
      out(idx_acq) = in;
    else
      out = in(idx_acq);
    end
  end

end


%-------------------------
%--- Support Functions ---
%-------------------------

function out = mri_reconSPIRiT_get_points( pt_idx, array, kernel_sz )
  % mri_reconSPIRiT_get_points helper function to retrieve correct points
  % This function is used both to set up the weights and get the appropriate points and
  % to interpolate with when filling in k-space.
  % It's slightly different than the GRAPPA version because of the points we need to grab
  %
  % What about if we know that the kernel is square and just need the size
  %
  % Author: Alex McManus
  % *********************
  %   Input Parameters:
  % *********************
  %
  %     pt_idx: A 3 element vector containing the [row column coilidx] of a specific
  %     point. Coil index matters here since for solving for the weights, we
  %     use the point at [row column] from the other coils.
  %         NOTE: this is different than grappa
  %
  %     array: the full 3D (ny x nx x ncoils) kspace data
  %
  %     kernel_sz: a 2 element array containing the [row column] size of the
  %     kernel. for spirit we'll always use this kernel completely filled in
  %     except for the center point
  %
  % *********************
  %   Output Variables:
  % ********************* 
  %
  %    out: the set of points for the input pt_idx over all the coils

  py = pt_idx(1);
  px = pt_idx(2);
  
  k_y = kernel_sz(1);
  k_x = kernel_sz(2);

  kdy = (k_y - 1) / 2;
  kdx = (k_x - 1) / 2;

  ypts = py-kdy:py+kdy;
  xpts = px-kdx:px+kdx;

  pts = array(ypts, xpts, :);
  out = pts(:);
end


function [W, A, B] = mri_reconSPIRiT_get_weights(acr, kernel_sz)
  %new_spirit_get_weights lets try this again
  
  % auto-calibration region dimensions
  szAcrY = size(acr, 1);
  szAcrX = size(acr, 2);
  ncoils = size(acr, 3);
  
  % kernel dimensions
  szKerY = kernel_sz(1);
  szKerX = kernel_sz(2);
  
  n = szKerY*szKerX;
  
  % centering
  center_y = (szKerY - 1)/2;
  center_x = (szKerX - 1)/2;
  
  % number of sliding kernel fits in the ACR
  numfits_x = szAcrX - szKerX + 1;
  numfits_y = szAcrY - szKerY + 1;
  
  npoints = szKerY*szKerX*ncoils;
  W = zeros([kernel_sz ncoils ncoils]);
  A = zeros([numfits_x*numfits_y npoints]);
  
  for coilIdx = 1:ncoils
    B = zeros( [numfits_x*numfits_y 1] );
    b_idx = 1;
    for x = 1:numfits_x
      for y = 1:numfits_y
        B(b_idx) = acr( y + center_y, x + center_x, coilIdx);
        pts = mri_reconSPIRiT_get_points([y+center_y x+center_x coilIdx], acr, kernel_sz);
        A(b_idx, :) = pts;
        b_idx = b_idx + 1;
      end
    end
  
    dk = zeros([szKerY szKerX ncoils]); 
    dk((end+1)/2, (end+1)/2, coilIdx) = 1;
    smp = ones(size(dk));
    smp( dk == 1 ) = 0;
    idxA = find(smp);
  
    A2 = A(:, idxA);
    rk = A2 \ B;
    
    ker = zeros(size(dk));
    ker(idxA) = rk;

    W(:, :, :, coilIdx) = ker;
  end
end


function out = spirit_conv(array, weights)
  nCoils = size( array, 3 );
  out = zeros( size( array ) );

  for coil = 1 : nCoils
    wi = weights( :, :, :, coil );
    for getCoil = 1 : nCoils
      wij = wi( :, :, getCoil );
      res = filter2( wij, array( :, :, getCoil ), 'same' );
      out( :, :, coil ) = out( :, :, coil ) + res;
    end
  end
end


function out = spirit_conv_adj( array, weights )
  nCoils = size( array, 3);
  out = zeros( size( array ) );

  for coil = 1 : nCoils
    wi = squeeze( weights( :, :, coil, : ) );
    for getCoil = 1 : nCoils
      wij = wi( :, :, getCoil );
      wij_conj = conj( wij );
      %ker = flipud( fliplr(wij_conj) );
      ker = rot90( wij_conj, 2 );
      res = filter2( ker, array( :, :, getCoil ), 'same' );
      out( :, :, coil ) = out( :, :, coil ) + res;
    end
  end
end

