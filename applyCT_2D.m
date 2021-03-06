
function out = applyCT_2D( f, kTraj, kCy, kCx, Cy, Cx, varargin )
  % out = applyCT_2D( f, kTraj, N, kCy, kCx, Cy, Cx [, newKs ] )
  %
  % Applies adjoint of a continuous circular convolution of a kernel as detailed in
  % http://nicholasdwork.com/tutorials/dworkGridding.pdf
  %
  % Inputs:
  % f - 
  % kTraj - 
  % kCy -
  % kCx -
  % Cy - 
  % Cx - 
  %
  % Optional Inputs:
  %  newKs - An nK x 2 array specifying the ky / kx (first / second column) coordinates
  %          of the new points
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
    newKs = size2fftCoordinates( [ size(f,1) size(f,2) ] );
    newKy = newKs{1};  newKx = newKs{2};
    gridKsSupplied = false;
  else
    newKy = newKs(:,1);  newKx = newKs(:,2);
    gridKsSupplied = true;
  end

  kTrajY = kTraj(:,1);    kyDistThresh = max( kCy );
  kTrajX = kTraj(:,2);    kxDistThresh = max( kCx );

  kyDists = min( abs( kTrajY - newKy' ), ...
                 abs( kTrajY + 1.0 - newKy' ) );
  kyDists = min( kyDists, ...
                 abs( kTrajY - 1.0 - newKy' ) );

  kxDists = min( abs( kTrajX - newKx' ), ...
                 abs( kTrajX + 1.0 - newKx' ) );
  kxDists = min( kxDists, ...
                 abs( kTrajX - 1.0 - newKx' ) );

	nTraj = size( kTraj, 1 );
  if gridKsSupplied == false

    nfs = size( f, 3 );
    segLength = 250;
    nSegs = ceil( nTraj / segLength );
    out = cell( nSegs, 1 );
    parfor segIndx = 1 : nSegs
      sIndx = ( segIndx - 1 ) * segLength + 1;
      eIndx = min( segIndx * segLength, nTraj );
      nSegTraj = eIndx - sIndx + 1;
      tmp = zeros( nSegTraj, nfs );

      for trajIndx = sIndx : eIndx
        shortIndxsY = find( kyDists( trajIndx, : ) < kyDistThresh );   %#ok<PFBNS>
        shortIndxsX = find( kxDists( trajIndx, : ) < kxDistThresh );   %#ok<PFBNS>

        CValsY = interp1( kCy, Cy, kyDists( trajIndx, shortIndxsY ), 'linear', 0 );
        CValsX = interp1( kCx, Cx, kxDists( trajIndx, shortIndxsX ), 'linear', 0 );
        CValsYX = CValsY' * CValsX;
        fCValsYX = bsxfun( @times, f( shortIndxsY, shortIndxsX, : ), CValsYX );
        tmp( trajIndx - sIndx + 1, : ) = reshape( sum( sum( fCValsYX, 1 ), 2 ), [ 1 nfs ] );
      end
      out{ segIndx } = tmp;
    end
    out = cell2mat( out );

  else

    shortIndxs = find( kyDists < kyDistThresh & kxDists < kxDistThresh );

    CValsY = interp1( kCy, Cy, kyDists( shortIndxs ), 'linear', 0 );
    CValsX = interp1( kCx, Cx, kxDists( shortIndxs ), 'linear', 0 );
    CValsYX = CValsY .* CValsX;

    nNew = size( newKs, 1 );
    [ trajShortIndxs, fShortIndxs ] = ind2sub( [ nTraj nNew ], shortIndxs );
    fThese = f( fShortIndxs, : );
    fCValsYX = bsxfun( @times, fThese, CValsYX );

    out = zeros( nTraj, size( f, 2 ) );
    for i = 1 : numel( shortIndxs )
      out( trajShortIndxs( i ), : ) = out( trajShortIndxs( i ), : ) + fCValsYX( i, : );
    end
  end

end
