
function figure = plotTriangles( triangles )  
  % figure = plotTriangles( triangles )
  %   creates shown 3D figure from given triangles
  % 
  % Inputs: 
  % triangles - 3D array of size 3 x 3 x nTriangles
  % 
  % Outputs: 
  % figure - 3D figure of polyhedron with triangle faces made from given
  %   triangle vertices
  %
  % Written by Jolie Wang - Copyright 2021
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  nTriangles = size( triangles, 3 );
  nPoints = size( triangles, 3 ) * 4;

  x = zeros( nPoints, 1 );
  y = zeros( nPoints, 1 );
  z = zeros( nPoints, 1 );

  index = 1;
  for i = 1:nTriangles

    xValues = triangles( :, 1, i );
    x( index ) = xValues( 1 );
    x( index + 1 ) = xValues( 2 );
    x( index + 2 ) = xValues( 3 );
    x( index + 3 ) = xValues( 1 );

    yValues = triangles( :, 2, i );
    y( index ) = yValues( 1 );
    y( index + 1 ) = yValues( 2 );
    y( index + 2 ) = yValues( 3 );
    y( index + 3 ) = yValues( 1 );

    zValues = triangles( :, 3, i );
    z( index ) = zValues( 1 );
    z( index + 1 ) = zValues( 2 );
    z( index + 2 ) = zValues( 3 );
    z( index + 3 ) = zValues( 1 );

    index = index + 4;
    
  end

  figure = plot3( x, y, z );
end
