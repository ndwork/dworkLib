function k_out = mri_reconGRAPPA( k_in, kernel_sz, acr_sz )
  % 
  % Generalized Autocalibrating Partially Parallel Acquisitions (GRAPPA) is a parallel imaging
  % reconstruction algorithm for magnetic resonance imaging. It reconstructs missing samples as a 
  % linear combination of nearby points in k-space from *all available coils* by using an 
  % autocalibration region (ACR), at the center of the image, to solve for the interpolation weights
  % given an interpolation kernel size.
  %
  % *********************
  %   Input Parameters:
  % *********************
  %
  %    k_in:  A 3D (size ny x nx x ncoils) array of measured k-space. The
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
  % *********************
  %   Output Variables:
  % *********************
  %
  %    k_out: A 3D (size ny x nx x ncoils) array of k-space values with the
  %    reconstructed values filled in.
  %
  % *********************
  %
  % Written by Alex McManus - Copyright 2022
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It is offered without any
  % warranty expressed or implied, including the implied warranties of merchantability or fitness
  % for a particular purpose.

  
  % kernel dimensions
  kernel_dy = kernel_sz(1);
  kernel_dx = kernel_sz(2);
  
  % auto-calibration region dimensions
  acr_dy = acr_sz(1);
  acr_dx = acr_sz(2);
  
  % size checking
  if mod(kernel_dx, 2) ~= 1 || mod(kernel_dy, 2) ~= 1
    error('Kernel size should be odd in both directions');
  end
  
  if mod(acr_dx, 2) ~= 1 || mod(acr_dy, 2) ~= 1
    error('ACR Size should be odd in both directions');
  end
  
  if kernel_dx > acr_dx
    error('Kernel is larger than ACR in kx-dimension');
  end
  
  if kernel_dy > acr_dy
    error('Kernel is larger than ACR in ky-dimension');
  end

  % get auto-calibration region
  acr = get_acr( k_in, acr_sz );

  % for a given size of kernel, we'll have edge cases depending on the
  % sampling pattern
  % for each of these edge cases, we'll have to solve for a different set of
  % weights
  % we assume the undersampling pattern is the same for each coil, so use the
  % first coil 
  k1 = squeeze(k_in(: ,:, 1));
  [kernels, karray] = get_kernels(k1, kernel_sz);
  k_out = k_in;

  for i = 1:numel(kernels)
    ka = kernels(i).ker;
    Wi = get_weights(acr, ka);
  
    karray_temp = (karray == i);
    oi = fill_points(k_in, karray_temp, ka, Wi);
    k_out = k_out + oi;
  end

end



% --------------------------------
%        Support functions
% --------------------------------

function out = fill_points(array, kernel_array, kernel, weights)
  %fill_points helper function to fill in missing k-space points
  % This function fills in the missing k-space data for a given kernel
  %
  % Author: Alex McManus
  % *********************
  %   Input Parameters:
  % *********************
  %
  %     array: the full 3D (ny x nx x ncoils) kspace data
  %
  %     kernel_array: a [ny x nx] array that is nonzero at points that use
  %     this interpolation kernel and zeros elsewhere
  %
  %     kernel: a 2D array of 1's and 0's that corresponds to the current
  %     GRAPPA kernel
  %
  %     weights: the previously solved-for interpolation weights for the
  %     given kernel
  %
  % *********************
  %   Output Variables:
  % ********************* 
  %
  %    out: a 3D (ny x nc x ncoils) array of reconstructed k-space data, only
  %    nonzero at the points specified in kernel_array
  
  
  [rows, cols] = find(kernel_array ~= 0);
  out = zeros(size(array));
  
  n = length(rows);
  for i = 1:n
    ri = rows(i);
    ci = cols(i);
    pts = get_points([ri ci], array, kernel);
  
    fillpt = weights.' * squeeze(pts);
    out(ri, ci, :) = fillpt;
  end

end


function acr = get_acr(k_in, acr_sz)
  % get_acr get ACR points from data
  % 
  % The point of this helper function is to retrieve the fully sampled region
  % in the middle of k-space - the auto-calibration region (ACR).
  % This takes the undersampled k-space data and the size of the ACR and
  % returns the values in the ACR to be used for solving for the weights in
  % GRAPPA or other algorithms
  %
  % Author: Alex McManus
  % *********************
  %   Input Parameters:
  % *********************
  %
  %    k_in:  A 3D (size ny x nx x ncoils) array of measured k-space. The
  %    assumption in GRAPPA is that there is a fully sampled region in the
  %    center of the image. The center point convention is as follows:
  %
  %    if our array is odd size, we choose the center point
  %           o o o x o o o
  %    if our array is even size, we choose length / 2 + 1
  %           o o o x o o 
  %
  %    acr_sz: A 2 element vector containing the dimensions of the ACR. ACR
  %    must be odd sized in both directions.
  %
  % *********************
  %   Output Variables:
  % *********************
  %
  %    acr: A 3D (size ky x kx x ncoils) array of k-space values for the
  %    auto-calibration region.
  
  ny = size(k_in, 1);
  nx = size(k_in, 2);
  
  acr_dy = acr_sz(1);
  acr_dx = acr_sz(2);
  
  if mod(acr_dx, 2) ~= 1 || mod(acr_dy, 2) ~= 1
    error('ACR Size should be odd in both directions');
  end
  
  acr_dx = (acr_dx - 1)/2;
  acr_dy = (acr_dy - 1)/2;
  
  % we need to make sure that we have the correct center
  % if our array is odd size, we choose the center point
  % o o o x o o o
  % if our array is even size, we choose length / 2 + 1
  % o o o x o o 
  center_y = ceil((ny+1)/2); 
  center_x = ceil((nx+1)/2);
  
  ypts = center_y - acr_dy : center_y + acr_dy;
  xpts = center_x - acr_dx : center_x + acr_dx;
  
  acr = k_in(ypts, xpts, :); 
  
  if sum(acr == 0, 'all') > 0
    error('Auto-Calibration Region not fully sampled.');
  end

end


function [out_struct, karray] = get_kernels(array, kernel_sz)
  % get_kernels find all of the kernels we'll need for a given kernel size
  %
  % We know in advance how large our kernel is. While most of the points
  % needing to be filled in will have the same set of points around them to
  % use, we will have edge cases (in most cases, literally at the edge of
  % collected k-space). This routine identifies all combinations of points
  % which will be used to fill in our reconstruction.
  %
  % Author: Alex McManus
  % *********************
  %   Input Parameters:
  % *********************
  %
  %    array:  A 2D (size ny x nx) array of measured k-space. The assumption
  %    here is that each coil will have the same undersampling pattern.
  %    Uncollected points are assumed to be zero.
  %
  %    kernel_sz: a 2 element vector containing the dimensions of the kernel.
  %    We use MATLAB dimensions, so it's [row column]
  %    For example, a 3x1 kernel will fit inside a 5x5 ACR as:
  %                 x o o o o
  %                 x o o o o
  %                 x o o o o
  %                 o o o o o 
  %                 o o o o o 
  %
  % *********************
  %   Output Variables:
  % *********************
  %
  %    out_struct: a structure that will have numel(out_struct) = # of unique
  %    kernels. out_struct(i) will be a matrix of dimensions kernel_sz that
  %    will have 1's at collected points and 0's at uncollected points.
  %
  %    This structure only has one field - ker. This is accessed as
  %           out_struct(i).ker
  %    which will be the aformentioned kernel matrix.
  %    
  %    karray: a 2D array of size (ny x nx). This will have 0's at
  %    *collected* points. A point with a nonzero value of i will use kernel
  %    out_struct(i).ker.
  
  karray = zeros(size(array));
  
  yrows = find(sum(array == 0, 2));
  xcols = find(sum(array == 0, 1));
  
  sKerY = kernel_sz(1);
  sKerX = kernel_sz(2);
  kcell = {};
  
  if mod(sKerY, 2) ~= 1 || mod(sKerX, 2) ~= 1
    error('Kernel dimensions must be odd');
  end
  
  ker_dy = (sKerY - 1) / 2;
  ker_dx = (sKerX - 1) / 2;
  
  for yi = 1:numel(yrows)
    y = yrows(yi);
    for xi = 1:numel(xcols)
      x = xcols(xi);
      if array(y, x) ~= 0
        continue
      end
      k1 = zeros([sKerY, sKerX]);
      for ky = -ker_dy:ker_dy
        if y + ky < 1 || y + ky > size(array, 1)
          continue
        end
  
        for kx = -ker_dx:ker_dx
          if x + kx < 1 || x + kx > size(array, 2)
            continue
          end
          if array(y+ky, x+kx) ~= 0
            k1(ky + ker_dy + 1, kx + ker_dx + 1) = 1;
          end
  
        end
      end
      
      n = numel(kcell);
      found = false;
      for cellidx = 1:n
        if kcell{cellidx} == k1
          found = true;
          karray(y, x) = cellidx;
          break;
        end
      end
  
      if ~found
        kcell(n+1) = {k1};
        karray(y, x) = n+1;
      end
  
    end
  end
  
  out_struct = struct( 'ker', kcell );
end


function out = get_points(pt_idx, array, kernel)
  % get_points helper function to retrieve correct points
  %
  % this function is used both to set up the weights and get the appropriate
  % points to interpolate with when filling in k-space
  %
  % Author: Alex McManus
  % *********************
  %   Input Parameters:
  % *********************
  %
  %     pt_idx: A 2 element vector containing the [row column] of a specific
  %     point. Coil index does not matter here since for solving for the
  %     interpolation weights and filling in the missing values, we need the
  %     values from *every* coil, so we pull from every coil all the time
  %     regardless of which coil we're *currently* solving for.
  %
  %     array: the full 3D (ny x nx x ncoils) kspace data
  %
  %     kernel: a 2D array of 1's and 0's that corresponds to the current
  %     GRAPPA kernel
  %
  % *********************
  %   Output Variables:
  % ********************* 
  %
  %    out: the set of points for the input pt_idx over all the coils
  
  k1 = logical(kernel);
  
  kdy = (size(kernel, 1) - 1) / 2;
  kdx = (size(kernel, 2) - 1) / 2;
  
  py = pt_idx(1);
  px = pt_idx(2);
  
  % need to add in error checking for the kernel here
  dist_ym = py - kdy;
  dist_yp = py + kdy;
  dist_xm = px - kdx;
  dist_xp = px + kdx;
  
  if dist_ym < 1
    % at the top of the array
    ypts = 1:py+kdy;
    cf = 1 - dist_ym; % correction factor
    k1 = k1(cf+1:end, :);
  elseif dist_yp > size(array, 1)
    % at the bottom
    ypts = py-kdy:size(array, 1);
    cf = dist_yp - size(array, 1);
    k1 = k1(1:end-cf, :);
  else
    ypts = py-kdy:py+kdy;
  end
  if dist_xm < 1
    % on the left edge
    xpts = 1:px+kdx;
    cf = 1 - dist_xm;
    k1 = k1(:, cf+1:end);
  elseif dist_xp > size(array, 2)
    % on the right edge
    xpts = px-kdx:size(array, 2);
    cf = dist_xp - size(array, 2);
    k1 = k1(:, 1:end-cf);
  else
    xpts = px-kdx:px+kdx;
  end
  
  pts = array(ypts, xpts, :);
  k2 = repmat(k1, [1, 1, size(array, 3)]);
  
  out = pts(k2);
end


function [W, A, B_total] = get_weights(acr, kernel)
  % get_weights solve for weights for given kernel
  %
  % This function solves for the prediction weights for a given kernel.
  % we'll solve B = A * W for W, where B is the collection of points of the
  % ACR, A are the corresponding points of the kernel, and W are the weights
  %
  % B = [ (numfits_x * numfits_y) x ncoils ]
  % W = [ (ncoils * npoints) x ncoils ]
  % A = [ (numfits_x * numfits_y) x (ncoils * npoints) ]
  %
  % Author: Alex McManus
  % *********************
  %   Input Parameters:
  % *********************
  %
  %    acr:  A 3D (size ny x nx x ncoils) array of fully sampled
  %    auto-calibration data.
  %
  %    kernel: a 2D array describing the current kernel. This is typically
  %    output from get_kernels.
  %
  %    For example, a 3x1 kernel will fit inside a 5x5 ACR as:
  %                 x o o o o
  %                 x o o o o
  %                 x o o o o
  %                 o o o o o 
  %                 o o o o o 
  %
  % *********************
  %   Output Variables:
  % *********************
  %
  %    W: a 2D array of dimensions [ (ncoils * npoints) ncoils ] containing
  %    the weights. The weights for coil i are W(:, i).
  %    
  %    (OPTIONAL) A: a 2D array of dimensions [ (numfits_x * numfits_y) x ncoils ] 
  %    containing the points used for each placement of the kernel within the
  %    ACR. Used for testing.
  %     
  %    (OPTIONAL) B: a 2D array of dimensions [ (numfits_x * numfits_y) (ncoils * npoints) ] 
  %    containing each of the center points for the kernel. Used for testing.
  
  % auto-calibration region dimensions
  acr_dy = size(acr, 1);
  acr_dx = size(acr, 2);
  ncoils = size(acr, 3);
  
  
  % kernel dimensions
  kernel_dy = size(kernel, 1);
  kernel_dx = size(kernel, 2);
  npoints = sum(kernel, 'all');
  
  % centering
  center_y = (kernel_dy - 1)/2;
  center_x = (kernel_dx - 1)/2;
  
  % number of sliding kernel fits in the ACR
  numfits_x = acr_dx - kernel_dx + 1;
  numfits_y = acr_dy - kernel_dy + 1;
  
  % we'll solve B = A * W for W, where B is the collection of points of the
  % ACR, A are the corresponding points of the kernel, and W are the weights
  
  % B = [ (numfits_x * numfits_y) x ncoils ]
  % W = [ (ncoils * npoints) x ncoils ]
  % A = [ (numfits_x * numfits_y) x (ncoils * npoints) ]
  
  % unroll the parallelization to just do one coil at a time
  W = zeros([ncoils*npoints ncoils]);
  B_total = zeros( [numfits_x * numfits_y 8]);
  for coilIndx = 1 : ncoils
    A = zeros( [numfits_x * numfits_y ncoils * npoints] );
    B = zeros( [numfits_x * numfits_y 1 ] );
    b_idx = 1;
    for kx = 1:numfits_x
      for ky = 1:numfits_y
        B(b_idx) = acr( ky+center_y, kx+center_x, coilIndx );
        ai = get_points([ky+center_y kx+center_x], acr, kernel);
        A(b_idx, :) = ai;
        b_idx = b_idx + 1;
      end
    end
  
    w = A \ B;
    W(:, coilIndx) = w;
    B_total(:, coilIndx) = B;
  end

end


