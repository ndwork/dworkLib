
function recon = mri_reconRoemer( coilRecons )
  % recon = mri_reconRoemer( coilRecons )
  %
  % Perform an optimal coil combination according to equation [32] of
  % "NMR Phased Array" by Roemer et al.
  %
  % Inputs:
  % coilRecons is an array of size ( Ny, Nx, nSlices, ..., nCoils ) of kSpace values
  %
  % Output:
  % recon is the reconstructed image
  %
  % Written by Nicholas Dwork - Copyright 2021
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp( 'Usage:  recon = mri_reconRoemer( kData [, ''multiSlice'', true/false ] )' );
    return;
  end

  ssqRecon = sqrt( sum( coilRecons .* conj( coilRecons ), ndims(coilRecons) ) );

  sMaps = bsxfun( @times, coilRecons, 1 ./ ssqRecon );

  recon = sum( bsxfun( @times, coilRecons, conj( sMaps ) ), ndims(coilRecons) );

end




