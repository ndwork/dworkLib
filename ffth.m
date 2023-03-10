
function out = ffth( in, varargin )
  % out = ffth( in )
  %
  % Compute the adjoint of the FFT
  %
  % Written by Nicholas Dwork, Copyright 2023
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addOptional( 'nFFT', [], @ispositive );
  p.addOptional( 'dim', [], @ispositive );
  p.parse( varargin{:} );
  nFFT = p.Results.nFFT;
  dim = p.Results.dim;

  if numel( dim ) > 0
    out = size( in, dim ) * ifft( in, nFFT, dim );
  else
    out = numel( in ) * ifft( in, nFFT, dim );
  end
end
