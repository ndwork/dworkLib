
function outPts = applyRadialDistortion2Pts( pts, ks, varargin )
  % out = applyRadialDistortion2Pts( pts, k[, c, 'dir', dir ] )
  %
  % Written according to section 7.4 of Multiple View Geometry, 2nd edition
  % by Hartley and Zisserman
  %
  % Inputs:
  % pts - an MxN array of points where M is the number of points and N
  %   is the dimension of the space
  % ks - the radial distortion coefficients
  %   L(r) = 1 + k(1) * r + k(2) * r^2 + ...
  % c - an N element array specifying the center of the image
  %
  % Optional Inputs:
  % dir - if 1, applies radial distortion.  If -1, undoes radial distortion
  % tol - tolerance to end fixed point iteration
  %
  % Outputs:
  % out - the image with radial distortion applied
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addOptional( 'c', zeros(2,1), @(x) isempty(x) || isnumeric(x) );
  p.addParameter( 'dir', 1, @(x) x == 1 || x == -1 );
  p.addParameter( 'maxIter', 50, @ispositive );
  p.addParameter( 'tol', 1d-10, @ispositive );
  p.parse( varargin{:} );
  c = p.Results.c;
  dir = p.Results.dir;
  maxIter = p.Results.maxIter;
  tol = p.Results.tol;

  rs = LpNorms( bsxfun( @minus, pts, c' ), 2, 2 );
  Ls = findLs( rs, ks );

  if dir == -1
    % outPts are the "corrected" undistorted points
    % Note: we don't want the Ls for the distorted values, which is what we would get by calling findLs here.
    % Instead, we want the Ls that exist for the undistorted values.  So we need to do an optimization.
    % Every radius relationship satisfies rUndistort = rDistort * L( rUndistort ).
    % Equivalently, rDistort = rUndistort / L( rUndistort ).
    % We solve for rUndistort using fixed point iteration

    rUndistort = rs; % initial guess
    for iter = 1 : maxIter
      betterLs = findLs( rUndistort, ks );
      betterLs( rs == 0 ) = 1;
      rUndistortNew = rs .* betterLs;

      dr = abs( rUndistortNew - rUndistort );
      rUndistort = rUndistortNew;
      if all( dr < tol ), break; end
    end

    outPts = bsxfun( @plus, bsxfun( @times, bsxfun( @minus, pts, c(:)' ), betterLs ), c(:)' );

  else
    % outPts are distorted with radial distortion
    outPts = bsxfun( @plus, bsxfun( @rdivide, bsxfun( @minus, pts, c(:)' ), Ls ), c(:)' );
  end
end


function Ls = findLs( rs, ks )
  rPower = rs;
  Ls = ones( size(rs) );
  for kIndx = 1 : numel( ks )
    Ls = Ls + ks(kIndx) .* rPower;
    if kIndx < numel( ks ), rPower = rPower .* rs; end
  end
end

