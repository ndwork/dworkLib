
function out = applyC_2D( F, kTraj, N, kCy, kCx, Cy, Cx, varargin )
  % out = applyC_2D( F, kTraj, N, kCy, kCx, Cy, Cx [, newKs ] )
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
  %  newKs - An nK x 2 array specifying the ky / kx (first / second column) coordinates
  %          of the new points
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
  p.addOptional( 'newKs', [] );
  p.parse( varargin{:} );
  newKs = p.Results.newKs;

  if numel( newKs ) == 0
    newKs = size2fftCoordinates( N );
    newKy = newKs{1};  newKx = newKs{2};
    gridKsSupplied = false;
  else
    newKy = newKs(:,1);  newKx = newKs(:,2);
    gridKsSupplied = true;
  end

  kTrajY = kTraj(:,1);
  kyDists = min( cat( 3, ...
    abs(   kTrajY         - newKy' ), ...
    abs( ( kTrajY + 1.0 ) - newKy' ), ...
    abs( ( kTrajY - 1.0 ) - newKy' ) ), [], 3 );
  CValsY = interp1( kCy, Cy, kyDists, 'linear', 0 );

  kTrajX = kTraj(:,2);
  kxDists = min( cat( 3, ...
    abs(   kTrajX         - newKx' ), ...
    abs( ( kTrajX + 1.0 ) - newKx' ), ...
    abs( ( kTrajX - 1.0 ) - newKx' ) ), [], 3 );
  CValsX = interp1( kCx, Cx, kxDists, 'linear', 0 );

  nFs = size( F,2 );
  if gridKsSupplied == false
    nTraj = size( kTraj, 1 );
    out = zeros( [ numel(newKy) numel(newKx) size(F,2) ] );
    for trajIndx = 1 : nTraj
      CValsYX = CValsY(trajIndx,:)' * CValsX(trajIndx,:);
      thisF = reshape( F(trajIndx,:), [ 1 1 nFs ] );
      out = out + bsxfun( @times, thisF, CValsYX );
    end

  else
    CValsYX = CValsY .* CValsX;
    out = zeros( numel(newKx), nFs );
    for indxF = 1 : nFs
      FCValsYX = bsxfun( @times, F(:,indxF), CValsYX );
      out(:,indxF) = transpose( sum( FCValsYX, 1 ) );
    end

  end

end
