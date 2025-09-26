
function [ intersectionPts, ts ] = lineCylinderIntersect( pt, vec, center, radius )
  % intersectionPts = lineCylinderIntersect( pt, vec, center, radius )
  %
  % Find the points of intersection between a line and a cylinder.  (Note that these points
  % may be imaginary or complex.)  The cylinder is assumed to be a circle when projected into the xy plane
  % with vertical extension in the z dimension
  %
  % The line is defined as { pt + t * vec, where t is real }
  % The cylinder is defined as
  %   { (x,y,z) : ( x - center(1) )^2 + ( y - center(2) )^2 = radius^2 }
  %
  % Outputs:
  % intersectionPts - an array specifying the intersection points
  % ts - an array of t values where the vec have norm 1
  %
  % Written by Nicholas Dwork, Copyright 2025
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 4
    disp( 'Usage: intersectionPts = lineCylinderIntersect( pt, vec, center, radius )' );
    if nargout > 0, intersectionPts = []; end
    if nargout > 1, ts = []; end
    return
  end

  if numel( center ) == 2, center = center(:)'; end
  if numel( pt ) == 2, pt = pt(:)'; end
  if numel( vec ) == 2, vec = vec(:)'; end

  if nargout > 1, vec = bsxfun( @rdivide, vec, LpNorms( vec, 2, 2 ) ); end

  [ ~, ts ] = lineCircleIntersect( pt(1:2), vec(:,[1 2]), center(1:2), radius );

  roots1Vec = bsxfun( @times, vec, ts(:,1) );
  roots2Vec = bsxfun( @times, vec, ts(:,2) );
  intersectionPts = [ bsxfun( @plus, transpose(pt), roots1Vec ); ...
                      bsxfun( @plus, transpose(pt), roots2Vec ); ];

end
