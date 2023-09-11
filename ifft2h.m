
function out = ifft2h( in )
  % out = ifft2h( in )
  %
  % Compute the adjoint of the two dimensional ifft
  %
  % Written by Nicholas Dwork, Copyright 2022
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if ismatrix( in )
    FHin = fft2( in );
  else
    FHin = fft( fft( in, [], 1 ), [], 2 );
  end

  out = FHin / ( size(in,1) * size(in,2) );
end
