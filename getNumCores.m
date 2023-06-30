
function [nPhysicalCores, nLogicalCores] = getNumCores()
  % [nPhysicalCores, nLogicalCores] = getNumCores()
  %
  % Gets the number of physical and logical cores available for processing with a parallel
  %   pool of workers (parpool).
  %
  % Outputs:
  % nPhysicalCores - the number of physical cores
  % nLogical Cores - the number of logical cores
  %
  % Written by Nicholas Dwork - Copyright 2023
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  coreInfo = evalc('feature(''numcores'')');
  newStr = splitlines(coreInfo);
  expression = '[0-9]+';

  expMatch = regexp( newStr{1}, expression, 'match' );
  nPhysicalCores = str2num( expMatch{1} );

  if nargout > 1
    expMatch = regexp( newStr{2}, expression, 'match' );
    nLogicalCores = str2num( expMatch{1} );
  end
end