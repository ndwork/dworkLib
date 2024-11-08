
function [out,sMaps] = mri_reconRoemer( in, varargin )
  % [out,sMaps] = mri_reconRoemer( in [, 'sMaps', sMaps, 'op', 'transp' or'notransp' ] )
  %
  % Perform an optimal coil combination according to equation [32] of
  %   "NMR Phased Array" by Roemer et al.
  %
  % Inputs:
  % in is an array of size ( Ny, Nx, nCoils, d1, d2, ... ) of kSpace values when op is 'notransp'
  %   Here, d1, ..., dN are optional dimensions
  %
  % Optional Inputs:
  % sMaps - array of sensitivity maps
  %
  % Output:
  % out is the reconstructed image when op is 'notransp'
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
    disp([ 'Usage:  recon = mri_reconRoemer( in [, ''sMaps'', sMaps, ', ...
           '''op'', ''transp'' or ''notransp'' ] )' ]);
    return;
  end

  p = inputParser;
  p.addParameter( 'op', [] );
  p.addParameter( 'sMaps', [], @isnumeric );
  p.parse( varargin{:} );
  op = p.Results.op;
  sMaps = p.Results.sMaps;

  if numel( op ) == 0, op = 'notransp'; end

  if numel( sMaps ) == 0
    if strcmp( op, 'transp' ), error( 'Must supply sMaps for transpose operation' ); end
    ssqRecon = sqrt( sum( in .* conj( in ), 3 ) );
    sMaps = bsxfun( @times, in, 1 ./ ssqRecon );
    sMaps( ~isfinite( sMaps ) ) = 0;
  end

  if strcmp( op, 'notransp' )
    out = sum( in .* conj( sMaps ), 3 );
  elseif strcmp( op, 'transp' )
    out = bsxfun( @times, sMaps, in );
  else
    error( 'Incorrect value of op supplied' );
  end
end

