
function out = uifft2( in, varargin )
  % out = uifft2( in, varargin )
  %
  % Compute the two-dimensional unitary inverse fft
  %
  % Optional Inputs:
  % Any optional inputs are passed to ifft2.
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

  out = sqrt(numel(in)) .* ifft2( in, varargin{:} );

end

