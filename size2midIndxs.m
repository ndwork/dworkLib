
function out = size2midIndxs( in )
  % out = size2midIndxs( in )
  %
  % Identify the middle indexes for a given size
  %
  % Inputs:
  % in - a 1D array specifying the size of the array in question in each dimension
  %
  % Outputs:
  % out - a 1D array specifying the middle indexes of the array in question
  %
  % Written by Nicholas Dwork - Copyright 2024
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  out = floor( 0.5 * in ) + 1;

end
