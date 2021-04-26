
function areas = calculateVoronoiAreas( kTraj, bounds )
  % areas = calculateVoronoiAreas( kTraj [, bounds ] )
  %
  % Calculate the areas of each Voronoi cell
  %
  % Optional Inputs:
  % bounds - [ xMin xMax yMin yMax ]

  if max( abs( imag( kTraj ) ) ) > 0
    kTraj = [ real( kTraj ) imag( kTraj ) ];
  end
  nTraj = size( kTraj, 1 );

  % uncomment these to plot voronoi diagram
  %[ vx, vy ] = voronoi( kTraj(:,2) , kTraj(:,1) );
  %plot( kTraj(:,2), kTraj(:,1), 'r.', vx, vy, 'b-' );
  %ylim([-0.6 0.6]);  xlim([-0.6 0.6]);  axis equal;

  % returns vertices and cells of voronoi diagram
  [ V, C ] = voronoin( kTraj ); 
  areas = zeros( nTraj, 1 );
  for trajIndx = 1 : numel( C )
    x = V( C{trajIndx}, 1 );  nx = numel( x );
    y = V( C{trajIndx}, 2 );  ny = numel( y );

    if nx == 0 || ny == 0
      A = NaN;

    elseif exist( 'bounds', 'var' ) && numel( bounds ) > 0 && ( ...
      min(x) < bounds(1) || max(x) > bounds(2) || ...
      min(y) < bounds(3) || max(y) > bounds(4) )
      A = NaN;

    else
      A = abs( sum( 0.5 * ( x( [2:nx 1] ) - x(:) ) .* ( y( [2:ny 1] ) + y(:) ) ) );

    end
    areas( trajIndx ) = A;
  end

end

