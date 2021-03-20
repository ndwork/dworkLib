function areas = voronoiAreas( pts )
  % areas = voronoiAreas( pts );
  %
  % input:  pts is a N x 2 array specifying the x and y locations of the points
  % output: area of cells for each point 
  %           (if point doesn't have neighbors the area is NaN)
  %
  % Written by Nicholas Dwork - Copyright 2019
  % Based on code written by John Pauly
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.


  % uncomment these to plot voronoi diagram
  %xs = pts(:,1);  ys = pts(:,2);
  %[vx, vy] = voronoi(xs,ys);
  %plot(xs,ys,'r.',vx,vy,'b-'); axis equal

  % returns vertices and cells of voronoi diagram
  nPts = size( pts, 1 );
  [ V, C ] = voronoin( pts ); 
  areas = zeros( nPts, 1 );
  for j = 1:length(pts)
    x = V(C{j},1);
    y = V(C{j},2);
    lxy = length(x);
    thisArea = abs(sum( 0.5*(x([2:lxy 1]) - x(:)).*(y([2:lxy 1]) + y(:))));
    areas(j) = thisArea;
  end


