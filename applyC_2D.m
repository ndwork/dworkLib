
function out = applyC_2D( F, traj, newTraj, kCy, kCx, Cy, Cx )
  % out = applyC_2D( F, traj, N, kCy, kCx, Cy, Cx )
  % or
  % out = applyC_2D( F, traj, newTraj, kCy, kCx, Cy, Cx )
  %
  % Applies a continuous circular convolution of a kernel as detailed in
  % http://nicholasdwork.com/tutorials/dworkGridding.pdf
  %
  % Inputs:
  %   F - 
  %   traj - 
  %   N - 
  %   newTraj - An nK x 2 array specifying the ky / kx (first / second column) coordinates
  %             of the new points
  %   kCy -
  %   kCx -
  %   Cy - 
  %   Cx - 
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if numel( newTraj ) == 2 && max( mod( newTraj, 1 ) ) == 0
    % newTraj is a two element array specifying the size of the grid
    N = newTraj;
    newTraj = size2fftCoordinates( N );
    newKy = newTraj{1};  newKx = newTraj{2};
    gridKsSupplied = true;
  else
    newKy = newTraj(:,1);  newKx = newTraj(:,2);
    gridKsSupplied = false;
  end

  kyDistThresh = max( kCy );
  kxDistThresh = max( kCx );

	nTraj = size( traj, 1 );
  if gridKsSupplied == true

    sOut = [ numel(newKy) numel(newKx) size(F,2) ];
    out = zeros( sOut );
    F = reshape( F, size(F,1), 1, size(F,2) );
    for trajIndx = 1 : nTraj

      kyDists = min( abs( traj(trajIndx,1)       - newKy ), ...
                     abs( traj(trajIndx,1) + 1.0 - newKy ) );
      kyDists = min( kyDists, ...
                     abs( traj(trajIndx,1) - 1.0 - newKy ) );

      kxDists = min( abs( traj(trajIndx,2)       - newKx ), ...
                     abs( traj(trajIndx,2) + 1.0 - newKx ) );
      kxDists = min( kxDists, ...
                     abs( traj(trajIndx,2) - 1.0 - newKx ) );

      shortIndxsY = find( kyDists < kyDistThresh );
      if numel( shortIndxsY ) == 0, continue; end

      shortIndxsX = find( kxDists < kxDistThresh );
      if numel( shortIndxsX ) == 0, continue; end

      CValsY = interp1( kCy, Cy, kyDists( shortIndxsY ), 'linear', 0 );
      CValsX = interp1( kCx, Cx, kxDists( shortIndxsX ), 'linear', 0 );
      CValsYX = CValsY * CValsX';

      outValues = bsxfun( @times, F(trajIndx,1,:), CValsYX );
      out( shortIndxsY, shortIndxsX, : ) = out( shortIndxsY, shortIndxsX, : ) + outValues;
    end

  else

    nNew = size( newTraj, 1 );
    nFs = size( F, 2 );
    
    segLength = 2000;
    nSegs = ceil( nTraj / segLength );

    out = cell( 1, 1, nSegs );
    parfor segIndx = 1 : nSegs
      startIndx = ( segIndx - 1 ) * segLength + 1;
      endIndx = min( segIndx * segLength, nTraj );

      tmp = zeros( nNew, nFs );
      for trajIndx = startIndx : endIndx
        kyDists = min( abs( traj(trajIndx,1)       - newKy ), ...
                       abs( traj(trajIndx,1) + 1.0 - newKy ) );   %#ok<PFBNS>
        kyDists = min( kyDists, ...
                       abs( traj(trajIndx,1) - 1.0 - newKy ) );

        kxDists = min( abs( traj(trajIndx,2)       - newKx ), ...
                       abs( traj(trajIndx,2) + 1.0 - newKx ) );
        kxDists = min( kxDists, ...
                       abs( traj(trajIndx,2) - 1.0 - newKx ) );

        shortIndxs = find( kyDists < kyDistThresh & kxDists < kxDistThresh );
        if numel( shortIndxs ) == 0, continue; end

        CValsY = interp1( kCy, Cy, kyDists( shortIndxs ), 'linear', 0 );
        CValsX = interp1( kCx, Cx, kxDists( shortIndxs ), 'linear', 0 );
        FCVals = bsxfun( @times, F( trajIndx, : ), CValsY .* CValsX );   %#ok<PFBNS>

        tmp( shortIndxs, : ) = tmp( shortIndxs, : ) + FCVals;
      end
      out{ segIndx } = tmp;
    end
    out = sum( cell2mat( out ), 3 );

  end

end
