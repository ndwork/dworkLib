
function out = ufft2( in )
  % out = ufft2( in )
  %
  % Compute the two-dimensional unitary fft of 
  % the first two dimensions of in
  %
  % Written by Nicholas Dwork, Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp( 'out = ufft2( in )' );
    if nargout > 0, out = []; end
    return
  end

  if ismatrix( in )
    Fin = fft2( in );
  else
    Fin = fft( fft( in, [], 1 ), [], 2 );
  end

  out = ( 1 / sqrt( size(in,1) * size(in,2) ) ) * Fin;
end
