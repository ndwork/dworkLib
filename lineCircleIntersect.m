
function [ intersectionPts, t ] = lineCircleIntersect( pt, vec, center, radius )
  % intersectionPts = lineCircleIntersect( pt, vec, center, radius )
  %
  % Find the points of intersection between a line and a cirlce.  (Note that these points
  % may be imaginary or complex.)
  % The line is defined as { pt + t * vec, where t is real }
  % The circle is defined as
  %   { (x,y) : ( x - center(:,1) )^2 + ( y - center(:,2) )^2 = radius^2 }
  %
  % Written by Nicholas Dwork, Copyright 2024
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if numel( pt ) == 2, pt = pt(:)'; end
  if numel( vec ) == 2, vec = vec(:)'; end
  if numel( center ) == 2, center = center(:)'; end

  a = vec(:,1) .* vec(:,1) + vec(:,2) .* vec(:,2);
  b = 2 * ( pt(:,1) .* vec(:,1) - vec(:,1) .* center(:,1) + pt(:,2) .* vec(:,2) - vec(:,2) .* center(:,2) );
  c = pt(:,1) .* pt(:,1) - 2 * pt(:,1) .* center(:,1) + center(:,1) .* center(:,1) ...
      + pt(:,2) .* pt(:,2) - 2 * pt(:,2) .* center(:,2) + center(:,2) .* center(:,2) - radius * radius;

  [roots1,roots2] = quadRoots( a, b, c );

  intersectionPts = [ pt + roots1(:) * vec; ...
                      pt + roots2(:) * vec; ];

  if nargout > 1, t = [ roots1(:); roots2(:); ]; end
end
