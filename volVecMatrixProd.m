
function out = volVecMatrixProd( vol, vec )
  % Each slice of vol is a matrix that will multiply a column vector
  %
  % Inputs:
  % vol - an MxNxK array where vol(:,:,k) left multiplies vec
  % vec - an Nx1 array that is left multiplied by each slice of vol
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

  out = vol(:,1,:) * vec(1);
  for i = 2:numel(vec)
    out = out + vec(i) * vol(:,i,:);
  end

end

