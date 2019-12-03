
function out = ufft2( in, varargin )
  % out = ufft2( in, varargin )
  %
  % Compute the two-dimensional unitary fft
  %
  % Optional Inputs:
  % Any optional inputs are passed to fft2.
  %
  % Written by Nicholas Dwork, Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if ~ismatrix( in ), error('Input must be two-dimensional' ); end

  out = 1/sqrt( size(in,1) * size(in,2) ) .* fft2( in, varargin{:} );

end
