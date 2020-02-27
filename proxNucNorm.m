
function out = proxNucNorm( in, thresh )
  % out = proxNucNorm( in, thresh )
  %
  % Returns the proximal operator of the nuclear norm of the input matrix
  %
  % Inputs:
  % in - an input matrix
  % thresh - the thresholding value
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if thresh == 0
    out = in;
    return;
  end

  [u,s,v] = svd( in, 'econ' );

  s = softThresh( diag(s), thresh );  % Singular values are always real

  out = u * diag(s) * v';

end
