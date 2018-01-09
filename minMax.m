
function out = minMax( in )
  % out = minMax( in )
  % Return the minimum and maximum of the input data as a two element array
  %
  % Written by Nicholas - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  out = [ min(in(:)) max(in(:)) ];
end