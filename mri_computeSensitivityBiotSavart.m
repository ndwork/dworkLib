
function sensitivities = mri_computeSensitivityBiotSavart( segs, locs )
  % sensitivities = mri_computeSensitivityBiotSavart( segs, locs )
  %
  % computes the sensitivity for each coil at all locations in locs
  % This code uses the technique described in "MRI image enhancement using 
  % Biot-Savart law at 3 Tesla" by Esin and Alpaslan
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
  bs = zeros( nLocs, 3 );

  s1 = segs( 1, : );
  for segIndx = 1 : nCoilSegs
    s2 = s1;
    s1 = segs( segIndx+1, : );
    s = s2 - s1;
    sUnit = s(:) / norm( s(:) );

    u = bsxfun( @plus, -locs, s2 );
    v = bsxfun( @minus, locs, s1 );

    cosAlphas = ( u * sUnit ) ./ norms( u, 2, 2 );
    cosBetas = ( v * sUnit ) ./ norms( v, 2, 2 );
    sinBetas = sin( acos( cosBetas ) );

    vNorms = norms( v, 2, 2 );
    Rs = vNorms ./ sinBetas;

    sMags = ( cosAlphas + cosBetas ) ./ Rs;  % sensitivity magnitudes

    sCross = makeCrossProdMatrix( sUnit );
    sDirs = ( sCross * v' )';
    sDirs = bsxfun( @times, sDirs, 1 ./ norms( sDirs, 2, 2 ) );

    bs = bs + bsxfun( @times, sDirs, sMags );
  end

  % By the principle of reciprocity, the sensitivity is proportional to the
  % magnitude of the emitted magnetic field at each location
  sensitivities = norms( bs, 2, 2 );
end
