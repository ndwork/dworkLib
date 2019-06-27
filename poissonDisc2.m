
function pts = poissonDisc2( r, varargin )
  % pts = poissonDisc2( r [, 'nK', nK, 'incSize', incSize, 'min_r', min_r ] )
  %
  % Written according to "Fast Poisson Disk Sampling in Arbitrary Dimensions" by Robert Bridson
  % and "Poisson Disk Sampling" written by Herman Tulleken located at
  % http://devmag.org.za/2009/05/03/poisson-disk-sampling/
  %
  % Inputs:
  % r - a scalar specifying the minimum size between points
  %     or a function handle for a function that accepts a point and returns a value
  %
  % Optional Inputs:
  % k - the number of points before rejection (default is 30)
  % incSize - vector memory allocation increment size (default is 5000)
  % min_r - if r is a function handle, then the user must also supply the minimum possible 
  %   r (must be positive)
  % verbose - if true, displays informative statements
  %
  % Outputs:
  % pts - a 2D array of size Num Points x 2
  %   each column represents the corresponding coordinate of the point
  %
  % Written by Nicholas Dwork, Copyright 2019
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addRequired( 'r' );
  p.addParameter( 'nK', 30, @ispositive );
  p.addParameter( 'incSize', 5000, @ispositive );
  p.addParameter( 'min_r', -1, @ispositive );
  p.addParameter( 'verbose', false, @(x) islogical(x) || isnumeric(x) );
  p.parse( r, varargin{:} );
  r = p.Results.r;
  nK = p.Results.nK;
  incSize = p.Results.incSize;
  min_r = p.Results.min_r;
  verbose = p.Results.verbose;

  b = [ [ -0.5, 0.5 ];     % vertical bounds
        [ -0.5, 0.5 ]; ];  % horizontal bounds

  nd = 2;  % Number of dimensions
  bDists = b(:,2) - b(:,1);

  pts = zeros( incSize, nd );
  pts_r = zeros( incSize, 1 );
  nPts = 1;
  pts(nPts,:) = rand(1,nd) .* bDists' + b(:,1)';
  if ispositive( r )
    pts_r(nPts) = r;
  else
    pts_r(nPts) = r( pts(1,:) );
  end

  if ispositive( r )
    min_r = min( r(:) );
  end
  cellSize = min_r / sqrt(nd);
  bGrid = cell( ceil( bDists' ./ cellSize ) );  % background grid

  gc = getGridCoordinate( pts(1,:), b, cellSize );
  bGrid{ gc(1), gc(2) } = 1;

  pList = zeros( incSize, 1 );  % processing list
  nProcPts = 1;
  pList( nProcPts ) = 1;
  
  while nProcPts > 0

    %-- Choose a point at random from the processing list
    pIndx = ceil( rand() * nProcPts );  % index into the processing list
    aIndx = pList( pIndx );  % active point index
    aPt = pts( aIndx, : );

    %-- Make k random points in the annulus between r and 2r from the active point
    active_r = pts_r( aIndx );
    kDists = active_r + active_r * rand(nK,1);
    kAngles = rand(nK,1) * 2*pi;
    kPts = [ aPt(1) + kDists .* cos( kAngles ),  ...
             aPt(2) + kDists .* sin( kAngles ) ];

    %-- Check to see if each candidate point is valid
    validPointFound = false;
    for kIndx = 1 : nK
      kPt = kPts( kIndx, : );
      if kPt(1) < b(1,1) || kPt(1) > b(1,2) || ...
         kPt(2) < b(2,1) || kPt(2) > b(2,2), continue; end

      if ispositive( r )
        this_r = r;
      else
        this_r = r( kPt );
      end

      % find the minimum distance to other points
      gc = getGridCoordinate( kPt, b, cellSize );
      nearbyIndxs = bGrid{ gc(1), gc(2) };
      nNearbyPts = numel( nearbyIndxs );
      if nNearbyPts > 0
        nearbyPts = pts( nearbyIndxs, : );
        tmp = nearbyPts - repmat( kPt, [ nNearbyPts 1 ] );
        dists2NearbyPts = sqrt( sum( tmp .* tmp, 2 ) );
        if min( dists2NearbyPts ) < this_r, continue; end
      end

      validPointFound = true;

      % add this point to the list of points
      nPts = nPts + 1;
      if size(pts,1) < nPts
        if verbose ~= false
          disp([ 'poissonDisc2: made ', num2str(nPts-1), ' points.' ]);
        end
        pts = [ pts; zeros(incSize,nd) ];   %#ok<AGROW>
        pts_r = [ pts_r; zeros(incSize,1) ];   %#ok<AGROW>
      end
      pts( nPts, : ) = kPt;
      pts_r( nPts ) = this_r;

      % add this point to the processing list
      nProcPts = nProcPts + 1;
      if numel( pList ) < nProcPts
        pList = [ pList; zeros(incSize,1); ];   %#ok<AGROW>
      end
      pList( nProcPts ) = nPts;

      % add this point to the grid in all squares within r distance
      gc = getGridCoordinate( kPt, b, cellSize );
      nNearCells = ceil( active_r / cellSize );
      for vIndx = max( gc(1) - nNearCells, 1 ) : min( gc(1) + nNearCells, size(bGrid,1) )
        for uIndx = max( gc(2) - nNearCells, 1 ) : min( gc(2) + nNearCells, size(bGrid,2) )
          bGrid{ vIndx, uIndx } = [ bGrid{ vIndx, uIndx }; nPts; ];
        end
      end
    end

    if validPointFound == false
      pList( pIndx ) = pList( nProcPts );
      nProcPts = nProcPts - 1;
    end
  end

  pts = pts( 1 : nPts, : );
end


function gc = getGridCoordinate( pt, bounds, cellSize )
  gc = [ ceil( ( pt(1) - bounds(1,1) ) / cellSize ), ...
         ceil( ( pt(2) - bounds(2,1) ) / cellSize ) ];
end

