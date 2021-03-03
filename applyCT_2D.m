
function out = applyCT_2D( f, kTraj, N, kCy, kCx, Cy, Cx, varargin )
  % out = applyCT_2D( f, traj, N, kCy, kCx, Cy, Cx [, gridKs, 'type', type ] )
  %
  % Applies the adjoint of a continuous circular convolution of a kernel as detailed in
  % http://nicholasdwork.com/tutorials/dworkGridding.pdf
  %
  % Inputs:
  % f - an array of size N
  % kTraj - Mx2 array where first/second column represents ky/kx location
  %         of trajectory
  % N - 2 element array specifying [Ny Nx] points in grid
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addOptional( 'gridKs', [] );
  p.parse( varargin{:} );
  gridKs = p.Results.gridKs;

  if numel( gridKs ) == 0
    gridKs = size2fftCoordinates( N );
    %[ gridKx, gridKy ] = meshgrid( gridKs{2}, gridKs{1} );
    gridKy = gridKs{1};
    gridKx = gridKs{2};
  else
    gridKy = gridKs(:,1);
    gridKx = gridKs(:,2);
  end

  kTrajY = kTraj(:,1);
  kyDists = min( cat( 3, ...
    abs(   kTrajY         - gridKy' ), ...
    abs( ( kTrajY + 1.0 ) - gridKy' ), ...
    abs( ( kTrajY - 1.0 ) - gridKy' ) ), [], 3 );
  CValsY = interp1( kCy, Cy, kyDists, 'linear', 0 );

  kTrajX = kTraj(:,2);
  kxDists = min( cat( 3, ...
    abs(   kTrajX         - gridKx' ), ...
    abs( ( kTrajX + 1.0 ) - gridKx' ), ...
    abs( ( kTrajX - 1.0 ) - gridKx' ) ), [], 3 );
  CValsX = interp1( kCx, Cx, kxDists, 'linear', 0 );

  nTraj = size( kTraj, 1 );
  out = zeros( nTraj, size(f,3) );
  for trajIndx = 1 : nTraj
    CValsYX = CValsY(trajIndx,:)' * CValsX(trajIndx,:);
    fCVals = bsxfun( @times, f, CValsYX );
    out(trajIndx,:) = sum( sum( fCVals, 2 ), 1 );
  end

  %nTraj = size( kTraj, 1 );
  %fCVals = zeros( nTraj, size(CValsY,2), size(CValsX,2) );
  %for trajIndx = 1 : nTraj
  %  fCVals(trajIndx,:,:) = f .* ( CValsY(trajIndx,:)' * CValsX(trajIndx,:) );
  %end
  %
  %out = squeeze( sum( sum( fCVals, 2 ), 3 ) );
end

