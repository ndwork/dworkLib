
function pts = poissonDisc2( r, varargin )
  % pts = poissonDisc2( r [, 'bIncSize', bIncSize, 'nK', nK, 'incSize', incSize, 'min_r', min_r ] )
  %
  % Written according to "Fast Poisson Disk Sampling in Arbitrary Dimensions" by Robert Bridson
  % and "Poisson Disk Sampling" written by Herman Tulleken located at
  % http://devmag.org.za/2009/05/03/poisson-disk-sampling/
  %
  % Inputs:
  % r - a scalar specifying the minimum size between points
  %     or a function handle for a function that accepts a point and returns a value
  %     or an image where each element of the image specifies the r value
  %
  % Optional Inputs:
  % bounds - a 2 x 2 array that specifies the bounds of the area where points get created
  %   The first/second row are the horizontal/vertical bounds, respectively.
  %   bounds(1,1) is the lower horizontal bound
  %   bounds(1,2) is the upper horizontal bound
  %   bounds(2,1) is the lower vertical bound
  %   bounds(2,2) is the upper vertical bound
  %   If r is an image and bounds is supplied then bounds must have the same aspect ratio as r
  % k - the number of points before rejection (default is 30)
  % incSize - vector memory allocation increment size (default is 5000)
  % min_r - if r is a function handle, then the user must also supply the minimum
  %   possible r (must be positive)
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
  p.addParameter( 'bounds', [], @isnumeric );
  p.addParameter( 'bIncSize', 20, @ispositive );
  p.addParameter( 'nK', 10, @ispositive );
  p.addParameter( 'incSize', 5000, @ispositive );
  p.addParameter( 'min_r', -1, @ispositive );
  p.addParameter( 'verbose', false, @(x) islogical(x) || isnumeric(x) );
  p.parse( r, varargin{:} );
  r = p.Results.r;
  bounds = p.Results.bounds;
  bIncSize = p.Results.bIncSize;
  nK = p.Results.nK;
  incSize = p.Results.incSize;
  min_r = p.Results.min_r;
  verbose = p.Results.verbose;

  nd = 2;  % Number of dimensions

  rIsFunction = true;
  if ispositive( r )
    min_r = min( r(:) );
    if numel( r ) > 1
      rImg = r;
      sImg = size( rImg );

      if numel( bounds ) > 0
        oldBounds = bounds;
        oldBDists = oldBounds(:,2) - oldBounds(:,1);
      end
      bounds = [ [ 1, sImg(2) ]; ...   % horizontal bounds
                 [ 1, sImg(1) ]; ];    % vertical bounds

      rScaling = mean( sImg(:) ./ oldBDists(:) );
      rImg = rImg * rScaling;

      min_r = min( rImg(:) );
      r = @(pts) rImg( round(pts(:,1)) + round( pts(:,2) - 1 ) * sImg(1) );
    else
      rIsFunction = false;
    end
  end

  if numel( bounds ) == 0
    bounds = [ [ -0.5, 0.5 ]; ...    % horizontal bounds
               [ -0.5, 0.5 ]; ];     % vertical bounds
  end
  bDists = bounds(:,2) - bounds(:,1);

  cellSize = min_r / sqrt(nd);
  bGrid = cell( ceil( bDists' ./ cellSize ) );  % background grid
  nBPts = zeros( ceil( bDists' ./ cellSize ) );  % number of points in each background grid cell
  tmp = zeros( bIncSize, 1 );
  for bIndx = 1 : numel( bGrid )
    bGrid{ bIndx } = tmp;
  end

  pts = zeros( incSize, nd );
  pts_r = zeros( incSize, 1 );
  nPts = 1;
  pts(nPts,:) = rand(1,nd) .* bDists' + bounds(:,1)';
  if rIsFunction
    pts_r(nPts) = r( pts(1,:) );
  else
    pts_r(nPts) = r;
  end

  gc = getGridCoordinate( pts(1,:), bounds, cellSize );
  nNearCells = ceil( pts_r(nPts) / cellSize );
  cL = gc(2) - nNearCells;  % left corner index
  cR = gc(2) + nNearCells;  % right corner index
  cB = gc(1) - nNearCells;  % bottom corner index
  cT = gc(1) + nNearCells;  % top corner index
  for uIndx = max( cL, 1 ) : min( cR, size(bGrid,2) )
    for vIndx = max( cB, 1 ) : min( cT, size(bGrid,1) )
      if ( uIndx == cL || uIndx == cR ) && ( vIndx == cB || vIndx == cT )
       continue
      end
      bGrid{ vIndx, uIndx }( 1 ) = nPts;
      nBPts( vIndx, uIndx ) = 1;
    end
  end
  maxNBPts = 1;  % Maximum number of points in any bGrid cell

  pList = zeros( incSize, 1 );  % processing list
  nProcPts = 1;
  pList( nProcPts ) = 1;

  if ~rIsFunction
    this_r = r;
    distThresh = this_r * this_r;
  end

  while nProcPts > 0

    %-- Choose a point at random from the processing list
    pIndx = ceil( rand() * nProcPts );  % index into the processing list
    aIndx = pList( pIndx );  % active point index
    aPt = pts( aIndx, : );
    active_r = pts_r( aIndx );

    validPointFound = false;

    %-- Make k random points in the annulus between r and 2r from the active point
    kDists = active_r + active_r * rand(nK,1);
    kAngles = rand(nK,1) * 2*pi;
    kPts = [ aPt(1) + kDists .* cos( kAngles ),  ...
            aPt(2) + kDists .* sin( kAngles ) ];
    kPts = kPts( kPts(:,1) > bounds(1,1) & kPts(:,1) < bounds(1,2) & ...
                 kPts(:,2) > bounds(2,1) & kPts(:,2) < bounds(2,2), : );

    if numel( kPts ) > 0
      kgcs = getGridCoordinate( kPts, bounds, cellSize );
      if rIsFunction
        krs = r( kPts );
      end

      %-- Check to see if each candidate point is valid
      for kIndx = 1 : size( kPts, 1 )
        kPt = kPts( kIndx, : );
        if rIsFunction
          this_r = krs( kIndx );
          distThresh = this_r * this_r;
        end

        % find the minimum distance to other points
        kgc = kgcs( kIndx, : );
        nNearbyPts = nBPts( kgc(1), kgc(2) );
        if nNearbyPts > 0
          nearbyIndxs = bGrid{ kgc(1), kgc(2) }( 1 : nNearbyPts );
          nearbyPts = pts( nearbyIndxs, : );
          nearbyPts(:,1) = nearbyPts(:,1) - kPt(1);
          nearbyPts(:,2) = nearbyPts(:,2) - kPt(2);
          dists2NearbyPtsSq = nearbyPts(:,1).^2 + nearbyPts(:,2).^2;
          if min( dists2NearbyPtsSq ) < distThresh, continue; end
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
        nNearCells = ceil( active_r / cellSize );
        cL = kgc(2) - nNearCells;  % left corner index
        cR = kgc(2) + nNearCells;  % right corner index
        cB = kgc(1) - nNearCells;  % bottom corner index
        cT = kgc(1) + nNearCells;  % top corner index

        for uIndx = max( cL, 1 ) : min( cR, size(bGrid,2) )
          for vIndx = max( cB, 1 ) : min( cT, size(bGrid,1) )          
            if ( uIndx == cL || uIndx == cR ) && ( vIndx == cB || vIndx == cT )
             continue
            end

            thisNBPts = nBPts( vIndx, uIndx );
            if maxNBPts >= bIncSize
              if thisNBPts == numel( bGrid{ vIndx, uIndx } )
                bGrid{ vIndx, uIndx } = [ bGrid{ vIndx, uIndx }; zeros( bIncSize, 1 ); ];
              end
            end
            thisNBPts = thisNBPts + 1;
            bGrid{ vIndx, uIndx }( thisNBPts ) = nPts;
            nBPts( vIndx, uIndx ) = thisNBPts;
            maxNBPts = max( maxNBPts, thisNBPts );
          end
        end
      end
    end

    if validPointFound == false
      pList( pIndx ) = pList( nProcPts );
      nProcPts = nProcPts - 1;
    end
  end

  if exist( 'oldBounds', 'var' )
    pts(:,1) = ( pts(:,1) - bounds(1,1) ) / rScaling + oldBounds(1,1);
    pts(:,2) = ( pts(:,2) - bounds(2,1) ) / rScaling + oldBounds(2,1);
  end

  pts = pts( 1 : nPts, : );
end


function out = ispositive( x )
  out = isnumeric(x) && min( x(:) > 0 );
end


function gc = getGridCoordinate( pts, bounds, cellSize )
  gc = ceil( [ ( pts(:,1) - bounds(1,1) ) , ...
               ( pts(:,2) - bounds(2,1) ) ] / cellSize );
end

