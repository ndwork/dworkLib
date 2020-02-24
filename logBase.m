
function out = logBase( in, base )
  % out = logBase( in, base )
  %
  % Compute the log to an arbitrary base
  %
  % Inputs:
  % in - an array of values
  % base - the base of the logarithm
  %
  % Written by Nicholas Dwork - Copyright 2020
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  out = log( in ) ./ log( base );

end

