
function out = applyC_2D( F, domTraj, rangeTraj, kCy, kCx, Cy, Cx )
  % out = applyC_2D( F, domTraj, N, kCy, kCx, Cy, Cx )
  % or
  % out = applyC_2D( F, N, rangeTraj, kCy, kCx, Cy, Cx )
  % or
  % out = applyC_2D( F, domTraj, rangeTraj, kCy, kCx, Cy, Cx )
  %
  % Applies a continuous circular convolution of a kernel as detailed in
  % http://nicholasdwork.com/tutorials/dworkGridding.pdf
  %
  % Inputs:
  %   F - An nTraj array representing the values of the function evaluated at each point in domTraj
  %   domTraj - An nTraj x 2 array specifying the ky / kx (first / second column) coordinates
  %             of the domain trajectory points
  %   rangeTraj - An nNew x 2 array specifying the ky / kx (first / second column) coordinates
  %             of the new points
  %             OR
  %             a two element array specifying the size of the grid in the Fourier domain
  %   kCy - array of convolution kernel domain values in y dimension
  %   kCx - array of convolution kernel domain values in x dimension
  %   Cy - array of convolution kernel values in y dimension
  %   Cx - array of convolution kernel values in x dimension
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if numel( domTraj ) == 2 && max( mod( domTraj, 1 ) ) == 0 && min( abs( domTraj ) ) > 0
    % domTraj is a two element array specifying the size of the grid
    N = domTraj;
    domTraj = size2fftCoordinates( N );
    domTrajKy = domTraj{1};  domTrajKx = domTraj{2};
    domTrajIsGrid = true;
  else
    domTrajKy = domTraj(:,1);  domTrajKx = domTraj(:,2);
    domTrajIsGrid = false;
  end

  if numel( rangeTraj ) == 2 && max( mod( rangeTraj, 1 ) ) == 0 && min( abs( rangeTraj ) ) > 0
    % rangeTraj is a two element array specifying the size of the grid
    N = rangeTraj;
    rangeTraj = size2fftCoordinates( N );
    rangeTrajKy = rangeTraj{1};  rangeTrajKx = rangeTraj{2};
    rangeTrajIsGrid = true;
  else
    rangeTrajKy = rangeTraj(:,1);  rangeTrajKx = rangeTraj(:,2);
    rangeTrajIsGrid = false;
  end

  if domTrajIsGrid == true && rangeTrajIsGrid == true
    error( 'This feature is not yet implemented.' );
  end

  kyDistThresh = max( kCy );
  kxDistThresh = max( kCx );

  segLength = 500;

	nTraj = size( domTraj, 1 );
  if rangeTrajIsGrid == true
    % rangeTraj is a grid and domTraj is not

    sF = size( F );
    nDomTraj = size( domTraj, 1 );
    sOut = [ numel( rangeTrajKy ) numel( rangeTrajKx ) prod( sF(2:end) ) ];
    F = reshape( F, [ nDomTraj prod( sF(2:end) ) ] );

    out = cell( numel( rangeTrajKy ), 1, 1 );

    parfor kyIndx = 1 : numel( rangeTrajKy )
      outSeg = zeros( 1, numel( rangeTrajKx ), prod( sF(2:end) ) );   %#ok<PFBNS>
      out{ kyIndx } = outSeg;

      kyDists = min( [ abs( rangeTrajKy( kyIndx ) - domTrajKy       ), ...
                       abs( rangeTrajKy( kyIndx ) - domTrajKy + 1.0 ), ...
                       abs( rangeTrajKy( kyIndx ) - domTrajKy - 1.0 ), ], [], 2 );

      shortIndxsY = find( kyDists < kyDistThresh );
      if numel( shortIndxsY ) == 0, continue; end

      CValsY = interp1( kCy, Cy, kyDists( shortIndxsY ), 'linear', 0 );

      subF = F( shortIndxsY, : );   %#ok<PFBNS>
      subDomTrajKx = domTrajKx( shortIndxsY );   %#ok<PFBNS>

      for kxIndx = 1 : numel( rangeTrajKx )
        kxDists = min( [ abs( rangeTrajKx( kxIndx ) - subDomTrajKx ), ...
                         abs( rangeTrajKx( kxIndx ) - subDomTrajKx + 1.0 ), ...
                         abs( rangeTrajKx( kxIndx ) - subDomTrajKx - 1.0 ), ], [], 2 );

        shortIndxsX = find( kxDists < kxDistThresh );
        if numel( shortIndxsX ) == 0, continue; end

        CValsX = interp1( kCx, Cx, kxDists( shortIndxsX ), 'linear', 0 );
        CValsYX = CValsY( shortIndxsX ) .* CValsX;

        outValues = bsxfun( @times, CValsYX, subF( shortIndxsX, : ) );
        outSeg( 1, kxIndx, : ) = sum( outValues, 1 );
      end
      out{ kyIndx } = outSeg;
    end

    out = cell2mat( out );
    out = reshape( out, sOut );

  elseif domTrajIsGrid == true
    % domTraj is a grid and rangeTraj is not

    sF = size( F );
    nF = prod( sF(3:end) );

    nRangeTraj = size( rangeTraj, 1 );
    F = reshape( F, [ sF(1:2) nF ] );

    nSegs = ceil( nRangeTraj / segLength );

    out = cell( nSegs, 1 );
    parfor segIndx = 1 : nSegs
      sIndx = ( segIndx - 1 ) * segLength + 1;
      eIndx = min( segIndx * segLength, nRangeTraj );
      nThisSeg = eIndx - sIndx + 1;
      segOut = zeros( nThisSeg, nF );

      for rTrajIndx = 1 : nThisSeg
        segOut( rTrajIndx, : ) = zeros( 1, nF );
        offsetIndx = ( segIndx - 1 ) * segLength;

        thisRangeTrajKy = rangeTrajKy( rTrajIndx + offsetIndx );   %#ok<PFBNS>
        kyDists = min( [ abs( thisRangeTrajKy - domTrajKy       ), ...
                         abs( thisRangeTrajKy - domTrajKy + 1.0 ), ...
                         abs( thisRangeTrajKy - domTrajKy - 1.0 ), ], [], 2 );

        shortIndxsY = find( kyDists < kyDistThresh );
        if numel( shortIndxsY ) == 0, continue; end

        thisRangeTrajKx = rangeTrajKx( rTrajIndx + offsetIndx );   %#ok<PFBNS>
        kxDists = min( [ abs( thisRangeTrajKx - domTrajKx ), ...
                         abs( thisRangeTrajKx - domTrajKx + 1.0 ), ...
                         abs( thisRangeTrajKx - domTrajKx - 1.0 ), ], [], 2 );

        shortIndxsX = find( kxDists < kxDistThresh );
        if numel( shortIndxsX ) == 0, continue; end

        CValsY = interp1( kCy, Cy, kyDists( shortIndxsY ), 'linear', 0 );
        CValsX = interp1( kCx, Cx, kxDists( shortIndxsX ), 'linear', 0 );
        CValsYX = CValsY * transpose( CValsX );

        subF = F( shortIndxsY, shortIndxsX, : );   %#ok<PFBNS>

        CF = bsxfun( @times, subF, CValsYX );
        segOut( rTrajIndx, : ) = reshape( sum( sum( CF, 1 ), 2 ), [ 1 nF ] );
      end
      out{ segIndx } = segOut;
    end
    out = cell2mat( out );

    if nF > 2
      out = reshape( out, [ nRangeTraj sF(3:end) ] );
    end

  else
    % Neither domTraj nor rangeTraj are a grid

    nRangeTraj = size( rangeTraj, 1 );
    nFs = size( F, 2 );

    nSegs = ceil( nTraj / segLength );

    out = cell( 1, 1, nSegs );
    parfor segIndx = 1 : nSegs
      startIndx = ( segIndx - 1 ) * segLength + 1;
      endIndx = min( segIndx * segLength, nTraj );

      tmp = zeros( nRangeTraj, nFs );
      for domTrajIndx = startIndx : endIndx
        kyDists = min( abs( domTraj(domTrajIndx,1)       - rangeTrajKy ), ...
                       abs( domTraj(domTrajIndx,1) + 1.0 - rangeTrajKy ) );   %#ok<PFBNS>
        kyDists = min( kyDists, ...
                       abs( domTraj(domTrajIndx,1) - 1.0 - rangeTrajKy ) );

        kxDists = min( abs( domTraj(domTrajIndx,2)       - rangeTrajKx ), ...
                       abs( domTraj(domTrajIndx,2) + 1.0 - rangeTrajKx ) );
        kxDists = min( kxDists, ...
                       abs( domTraj(domTrajIndx,2) - 1.0 - rangeTrajKx ) );

        shortIndxs = find( kyDists < kyDistThresh & kxDists < kxDistThresh );
        if numel( shortIndxs ) == 0, continue; end

        CValsY = interp1( kCy, Cy, kyDists( shortIndxs ), 'linear', 0 );
        CValsX = interp1( kCx, Cx, kxDists( shortIndxs ), 'linear', 0 );
        FCVals = bsxfun( @times, F( domTrajIndx, : ), CValsY .* CValsX );   %#ok<PFBNS>

        tmp( shortIndxs, : ) = tmp( shortIndxs, : ) + FCVals;
      end
      out{ segIndx } = tmp;
    end
    out = sum( cell2mat( out ), 3 );

  end

end
