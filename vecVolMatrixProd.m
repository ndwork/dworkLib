
function out = vecVolMatrixProd( vec, vol )
  % Each slice of vol is a matrix that will be left multiplied a row vector
  %
  % Inputs:
  % vec - an 1xM array that left multiplies each slice of vol
  % vol - an MxNxK array where vol(:,:,k) left multiplies vec
  %
  % Outputs:
  % out - an Mx1xK array
  %
  % Written by Nicholas Dwork, Copyright 2018
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  out = vec(1) * vol(1,:,:);
  for i = 2:numel(vec)
    out = out + vec(i) * vol(i,:,:);
  end

end

