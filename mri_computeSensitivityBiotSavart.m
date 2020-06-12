
function sensitivities = mri_computeSensitivityBiotSavart( segs, locs )
  % sensitivities = mri_computeSensitivityBiotSavart( segs, locs )
  %
  % computes the sensitivity for each coil at all locations in locs
  % This code uses the technique described in "MRI image enhancement using 
  % Biot?Savart law at 3 tesla" by Esin and Alpaslan
  %
  % Inputs:
  % segs - The coil is specified as a 2D array that correspond to a set of 3D
  %   coordinates (x,y,z).  The array is size Nx3, where N is the number of pieces
  %   of metal in the coils.  It is assumed that the coil is made of N-1 straight, 
  %   and that each piece are connected by coordinates one row in the array and its
  %   subsequent row.  If the coil is to be closed, then the last row and the first
  %   row should be the same coordinate.
  %   If the coil locations are an Nx2 array, it is assumed that the coil is located
  %   in the z=0 plane.
  %
  % locs - a 2D array of size Mx3 specifying the locations wherewe would like to
  %   compute the sensitivity.
  %
  % Outputs:
  % sensitivities - a 1D array that corresponds to the sensitivity at each location
  %   in locs
  %
  % Written by Nicholas Dwork, Copyright 2020
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  sCoil = size( segs );
  nCoilSegs = sCoil(1)-1;
  if nCoilSegs < 1, error( 'You did not specify enough coil segments' ); end
  if sCoil(2) == 2, segs = [ segs zeros( nCoilSegs ) ]; end
  
  nLocs = size( locs, 1 );
  sensitivities = zeros( nLocs, 1 );

  s1 = segs( 1, : );
  for segIndx = 1 : nCoilSegs
    s2 = s1;
    s1 = segs( segIndx, : );
    s = s2 - s1;
    sUnit = s / norm( s(:) );

    for locIndx = 1 : nLocs
      p = locs( locIndx, : );
      u = s2 - p;
      v = p - s1;

      cosAlpha = dotP( sUnit, u ) / norm( u(:) );
      cosBeta = dotP( sUnit, v ) / norm( v(:) );
      R = p - projV1ontoV2( v, s );

      sensitivities( locIndx ) = ( cosAlpha + cosBeta ) / R;
    end

  end


end

