
function out = applyC_2D( F, kTraj, N, kCy, kCx, Cy, Cx, varargin )
  % out = applyC_2D( F, kTraj, N, kCy, kCx, Cy, Cx [, newKs ] )
  %
  % Applies a continuous circular convolution of a kernel as detailed in
  % http://nicholasdwork.com/tutorials/dworkGridding.pdf
  %
  % Inputs:
  % F - 
  % kTraj - 
  % N - 
  % kCy -
  % kCx -
  % Cy - 
  % Cx - 
  %
  % Optional Inputs:
  %  newKs - An nK x 2 array specifying the ky / kx (first / second column) coordinates
  %          of the new points
  %  type -  by default, performs a circular convolution  If type=='noCirc',
  %          then it performs a regular(non-circular) convolution.
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addOptional( 'newKs', [] );
  p.parse( varargin{:} );
  newKs = p.Results.newKs;

  if numel( newKs ) == 0
    newKs = size2fftCoordinates( N );
    newKy = newKs{1};  newKx = newKs{2};
    gridKsSupplied = false;
  else
    newKy = newKs(:,1);  newKx = newKs(:,2);
    gridKsSupplied = true;
  end

  kTrajY = kTraj(:,1);    kyDistThresh = max( kCy );
  kTrajX = kTraj(:,2);    kxDistThresh = max( kCx );

  kyDists = min( abs( kTrajY - newKy' ), ...
                 abs( kTrajY + 1.0 - newKy' ) );
  kyDists = min( kyDists, ...
                 abs( kTrajY - 1.0 - newKy' ) );

  kxDists = min( abs( kTrajX - newKx' ), ...
                 abs( kTrajX + 1.0 - newKx' ) );
  kxDists = min( kxDists, ...
                 abs( kTrajX - 1.0 - newKx' ) );

	nTraj = size( kTraj, 1 );
  if gridKsSupplied == false

    sOut = [ numel(newKy) numel(newKx) size(F,2) ];
    out = zeros( sOut );
    F = reshape( F, size(F,1), 1, size(F,2) );
    for trajIndx = 1 : nTraj
      
      shortIndxsY = find( kyDists(trajIndx,:) < kyDistThresh );
      shortIndxsX = find( kxDists(trajIndx,:) < kxDistThresh );

      CValsY = interp1( kCy, Cy, kyDists( trajIndx, shortIndxsY ), 'linear', 0 );
      CValsX = interp1( kCx, Cx, kxDists( trajIndx, shortIndxsX ), 'linear', 0 );
      CValsYX = CValsY' * CValsX;

      outValues = bsxfun( @times, F(trajIndx,1,:), CValsYX );
      out( shortIndxsY, shortIndxsX, : ) = out( shortIndxsY, shortIndxsX, : ) + outValues;
    end

  else

    shortIndxs = find( kyDists < kyDistThresh & kxDists < kxDistThresh );

    CValsY = interp1( kCy, Cy, kyDists( shortIndxs ), 'linear', 0 );
    CValsX = interp1( kCx, Cx, kxDists( shortIndxs ), 'linear', 0 );
    CValsYX = CValsY .* CValsX;
    FCValsYX = bsxfun( @times, F( mod( shortIndxs-1, nTraj ) + 1, : ), CValsYX );

    out = zeros( size( kxDists ) );
    out( shortIndxs ) = FCValsYX;
    out = sum( out, 2 );
  end

end
