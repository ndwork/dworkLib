
function out = makeC_2D( domTraj, rangeTraj, kCy, kCx, Cy, Cx )
  % out = makeC_2D( domTraj, N, kCy, kCx, Cy, Cx )
  % or
  % out = makeC_2D( N, rangeTraj, kCy, kCx, Cy, Cx )
  % or
  % out = makeC_2D( domTraj, rangeTraj, kCy, kCx, Cy, Cx )
  %
  % Makes the circular convolution matrix C as detailed in
  % http://nicholasdwork.com/tutorials/dworkGridding.pdf
  %
  % Inputs:
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
  % Written by Nicholas Dwork - Copyright 2022
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

    error( 'Not yet implemented' );

  elseif domTrajIsGrid == true
    % domTraj is a grid and rangeTraj is not

    error( 'Not yet implemented' );

  else
    % Neither domTraj nor rangeTraj are a grid
    nSegs = ceil( nTraj / segLength );

    rows = cell( nSegs, 1 );
    cols = cell( nSegs, 1 );
    vals = cell( nSegs, 1 );

    parfor segIndx = 1 : nSegs
      startIndx = ( segIndx - 1 ) * segLength + 1;
      endIndx = min( segIndx * segLength, nTraj );

      nIncrement = 10000;
      nSegVals = nIncrement;
      segRows = zeros( nSegVals, 1 );
      segCols = zeros( nSegVals, 1 );
      segVals = zeros( nSegVals, 1 );
      nVals = 0;

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
        nShortDists = numel( shortIndxs );
        if nShortDists == 0, continue; end

        CValsY = interp1( kCy, Cy, kyDists( shortIndxs ), 'linear', 0 );
        CValsX = interp1( kCx, Cx, kxDists( shortIndxs ), 'linear', 0 );

        if nVals + nShortDists > nSegVals
          tmpVals = zeros( nSegVals + nShortDists + nIncrement, 1 );
          tmpVals( 1 : nVals ) = segVals( 1 : nVals );
          segVals = tmpVals;
        end

        segRows( nVals + 1 : nVals + nShortDists ) = domTrajIndx;
        segCols( nVals + 1 : nVals + nShortDists ) = shortIndxs;
        segVals( nVals + 1 : nVals + nShortDists ) = CValsY .* CValsX;
        nVals = nVals + nShortDists;
      end

      rows{ segIndx } = segRows( 1 : nVals );
      cols{ segIndx } = segCols( 1 : nVals );
      vals{ segIndx } = segVals( 1 : nVals );
    end

    rows = cell2mat( rows );
    cols = cell2mat( cols );
    vals = cell2mat( vals );
  end

  out = sparse( rows, cols, vals );
end
