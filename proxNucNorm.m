
function out = proxNucNorm( in, thresh )
  % out = proxL1Complex( in, thresh )
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

  [u,s,v] = svd( in );
  
  if isreal( in )
    s = softThresh( diag(s), thresh );
  else
    s = proxL1Complex( diag(s), thresh );
  end

  out = u * diag(s) * v;
end
