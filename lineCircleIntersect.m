
function [ intersectionPts, ts ] = lineCircleIntersect( pt, vec, center, radius )
  % intersectionPts = lineCircleIntersect( pt, vec, center, radius )
  %
  % Find the points of intersection between a line and a cirlce.  (Note that these points
  % may be imaginary or complex.)
  % The line is defined as { pt + t * vec, where t is real }
  % The circle is defined as
  %   { (x,y) : ( x - center )^2 + ( y - center )^2 = radius^2 }
  %
  % Outputs:
  % intersectionPts - an array specifying the intersection points
  % ts - an array of t values where the vec have norm 1
  %
  % Written by Nicholas Dwork, Copyright 2024
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 4
    disp( 'Usage: intersectionPts = lineCircleIntersect( pt, vec, center, radius )' );
    if nargout > 0, intersectionPts = []; end
    if nargout > 1, ts = []; end
    return
  end

  if numel( center ) == 2, center = center(:)'; end
  if numel( pt ) == 2, pt = pt(:)'; end
  if numel( vec ) == 2, vec = vec(:)'; end

  if nargout > 1, vec = bsxfun( @rdivide, vec, LpNorms( vec, 2, 2 ) ); end

  a = vec(:,1) .* vec(:,1) + vec(:,2) .* vec(:,2);

  if min( size( pt ) == size( vec ) ) == 1 && min( size( pt ) == size( center ) )
    b = 2 * ( pt(:,1) .* vec(:,1) - vec(:,1) .* center(:,1) + pt(:,2) .* vec(:,2) - vec(:,2) .* center(:,2) );
    c = pt(:,1) .* pt(:,1) - 2 * pt(:,1) .* center(:,1) + center(:,1) .* center(:,1) ...
        + pt(:,2) .* pt(:,2) - 2 * pt(:,2) .* center(:,2) + center(:,2) .* center(:,2) - radius(:) .* radius(:);

  else
    if min( size( pt ) == size( vec ) ) == 1
      ptVec11 = pt(:,1) .* vec(:,1);
      ptVec22 = pt(:,2) .* vec(:,2);
    else
      if numel( pt ) == 2
        ptVec11 = pt(1) * vec(:,1);
        ptVec22 = pt(2) * vec(:,2);
      else
        ptVec11 = pt(:,1) * vec(1);
        ptVec22 = pt(:,2) * vec(2);
      end
    end

    if min( size( vec ) == size( center ) ) == 1
      if numel( vec ) == 2, vec = vec(:)'; end

      vecCenter11 = vec(:,1) .* center(:,1);
      vecCenter22 = vec(:,2) .* center(:,2);
    else
      if numel( vec ) == 2
        vecCenter11 = vec(1) * center(:,1);
        vecCenter22 = vec(2) * center(:,2);
      else
        vecCenter11 = vec(:,1) * center(1);
        vecCenter22 = vec(:,2) * center(2);
      end
    end

    b = 2 * ( ptVec11 - vecCenter11 + ptVec22 - vecCenter22 );

    if min( size( pt ) == size( center ) ) == 1
      if numel( pt ) == 2, pt = pt(:)'; end
      if numel( center ) == 2, center = center(:)'; end

      ptCenter11 = pt(:,1) .* center(:,1);
      ptCenter22 = pt(:,2) .* center(:,2);
    else
      if numel( pt ) == 2
        ptCenter11 = pt(1) * center(:,1);
        ptCenter22 = pt(2) * center(:,2);
      else
        ptCenter11 = pt(:,1) * center(1);
        ptCenter22 = pt(:,2) * center(2);
      end
    end

    ptPt11 = pt(:,1) .* pt(:,1);
    ptPt22 = pt(:,2) .* pt(:,2);
    centerCenter11 = center(:,1) .* center(:,1);
    centerCenter22 = center(:,2) .* center(:,2);

    cPart1 = bsxfun( @minus, ptPt11, 2 * ptCenter11 );
    cPart1 = bsxfun( @plus, cPart1, centerCenter11 );
    cPart2 = bsxfun( @minus, ptPt22, 2 * ptCenter22 );
    cPart2 = bsxfun( @plus, cPart2, centerCenter22 );
    c = bsxfun( @minus, cPart1 + cPart2, radius(:) .* radius(:) );

  end

  [roots1,roots2] = quadRoots( a, b, c );

  roots1Vec = bsxfun( @times, roots1, vec );
  roots2Vec = bsxfun( @times, roots2, vec );
  intersectionPts = [ bsxfun( @plus, pt, roots1Vec ); ...
                      bsxfun( @plus, pt, roots2Vec ); ];

  if nargout > 1, ts = [ roots1(:); roots2(:); ]; end
end
