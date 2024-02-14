
function [recon,sMaps] = mri_reconRoemer( coilRecons, varargin )
  % [recon,sMaps] = mri_reconRoemer( coilRecons [, 'sMaps', sMaps ] )
  %
  % Perform an optimal coil combination according to equation [32] of
  %   "NMR Phased Array" by Roemer et al.
  %
  % Inputs:
  % coilRecons is an array of size ( Ny, Nx, ..., nCoils ) of kSpace values
  %
  % Optional Inputs:
  % sMaps - array of sensitivity maps
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
    disp( 'Usage:  recon = mri_reconRoemer( coilRecons [, ''sMaps'', sMaps ] )' );
    return;
  end

  p = inputParser;
  p.addParameter( 'sMaps', [], @isnumeric );
  p.parse( varargin{:} );
  sMaps = p.Results.sMaps;

  if numel( sMaps ) == 0
    ssqRecon = sqrt( sum( coilRecons .* conj( coilRecons ), ndims(coilRecons) ) );
    sMaps = bsxfun( @times, coilRecons, 1 ./ ssqRecon );
    sMaps( ~isfinite( sMaps ) ) = 0;
  end

  recon = sum( bsxfun( @times, coilRecons, conj( sMaps ) ), ndims(coilRecons) );
end

