
function out = fftshift2( in )
  % out = fftshift2( in )
  %
  % Performs fftshift on the first two dimensions of in
  %
  % Inputs:
  % in - input array of at least two dimensions
  %
  % Outputs:
  % out - output array
  %
  % Written by Nicholas Dwork - Copyright 2021
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  if nargin < 1
    disp( 'Usage:  out = fftshift2( in )' );
    if nargout > 0, out = []; end
    return;
  end

  out = fftshift( fftshift( in, 1 ), 2 );
end
