
function out = nSigDims( arrayIn )
  % out = nSigDims( arrayIn )
  % Determines the number of dimensions with size greater than 1
  %
  % Input:
  % arrayIn - the array to be evaluated
  %
  % Output:
  % out - the number of significant dimensions
  %
  % Written by Nicholas Dwork - Copyright 2017
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  s = size( arrayIn );
  out = sum( s > 1 );
end