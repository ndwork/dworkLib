
function out = applyC_2D( F, kTraj, N, kCy, kCx, Cy, Cx, varargin )
  % out = applyC_2D( F, kTraj, N, kCy, kCx, Cy, Cx [, gridKs, 'type', type ] )
  %
  % Applies a continuous circular convolution of a kernel as detailed in
  % http://nicholasdwork.com/tutorials/dworkGridding.pdf
  %
  % Inputs:
  % F - 
  % kTraj - 
  % N - 
  % kCy -
  % kCx -
  % Cy - 
  % Cx - 
  %
  % Optional Inputs:
  %  gridKs - 
  %  type -  by default, performs a circular convolution  If type=='noCirc',
  %          then it performs a regular(non-circular) convolution.
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
    gridKy = gridKs{1};  gridKx = gridKs{2};
  else
    gridKy = gridKs(:,1);  gridKx = gridKs(:,2);
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
  out = zeros( [ N size(F,2) ] );
  for trajIndx = 1 : nTraj
    CValsYX = CValsY(trajIndx,:)' * CValsX(trajIndx,:);
    thisF = reshape( F(trajIndx,:), [ 1 1 size(F,2) ] );
    out = out + bsxfun( @times, thisF, CValsYX );
  end

end
