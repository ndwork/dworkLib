function out = indicatorFunction( in, S )
  % out = indicatorFunction( in, S )
  %
  % outputs 0 if all the elments of the input are within S and Inf otherwise
  %
  % Inputs:
  % in - an array (of any size)
  %
  % Optional Inputs:
  % S - A two element array specifying the interval for the indicator function
  %     If a single sided bound is desired, set S with Inf.  For example, if
  %     in must be non-negative, set S to [ 0, Inf ]
  %     default is [ 0 0 ]
  %
  % Outputs:
  % out - 0 or Inf
  %
  % Written by Nicholas Dwork - Copyright 2021
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 2, S = [0 0]; end

  if min( in ) >= min( S ) && max( in ) <= max( S )
    out = 0;
  else
    out = Inf;
  end
end
