
function cellStrs = nums2cellStrs( numArray )
  % cellStrs = nums2cellStrs( numArray )
  %
  % Converts an array of numeric values to a cell array of strings
  %
  % Written by Nicholas - Copyright 2018
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp('Usage: cellStrs = nums2cellStrs( numArray )');
    return
  end

  cellNums = num2cell( numArray );
  cellStrs = cellfun( @num2str, cellNums, 'uni', false );

end