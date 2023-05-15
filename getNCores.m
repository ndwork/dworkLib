
function nCores = getNCores

  % nCores = getNCores()
  %
  % Get the number of logical cores available on the computer for parallel processing
  %
  % Written by Nicholas Dwork - Copyright 2023
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  nCores = [];
  evalc('nCores = feature(''numcores'')');
end