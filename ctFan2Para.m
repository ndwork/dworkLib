
function out = ctFan2Para( sinogram, gammas, betas, d, oThetas, oLines )
  % out = ctFan2Para( sinogram, gammas, betas, d, oThetas, oLines )
  %
  % Interpolates a sinogram created by a fan projection into a parallel projection
  % sinogram using physical units of a computed tomography detector
  %
  % Inputs:
  % gammas - a two element array specifying the min and max angle of the fan beam
  % betas - a two element array specifying the min and max rotation angles
  % d - the distance between the source and the isocenter
  % oBetas - a three element array specifying the min, delta, and max output rotation angles
  % oLines - the min, delta, and max spacing of projection lines
  %
  % Written by Nicholas Dwork - Copyright 2013
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  nThetas = round( ( oThetas(3) - oThetas(1) ) / oThetas(2) );
  nLines  = round( (  oLines(3) -  oLines(1) ) /  oLines(2) );

  thetas1 = linspace( oThetas(1), oThetas(3), nThetas );
  R1 = linspace( oLines(1), oLines(3), nLines );

  thetas = thetas1' * ones(1,nLines);
  Rs = ones(nThetas,1) * R1;

  oGammas = asin( Rs ./ d );
  oBetas = thetas + oGammas - pi/2;

  minBetas = min( betas(:) );
  maxBetas = max( betas(:) );
  indxs = find( oBetas > maxBetas );
  if numel( indxs ) > 0, oBetas( indxs ) = oBetas( indxs ) - 2*pi; end
  indxs = find( oBetas < minBetas );
  if numel( indxs ) > 0, oBetas( indxs ) = oBetas( indxs ) + 2*pi; end
  oBetas = min( oBetas, maxBetas );

  out = interp2( gammas, betas, sinogram, oGammas, oBetas, 'linear', 0 );
end
