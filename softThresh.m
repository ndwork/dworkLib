
function out = softThresh( in, thresh )
  % out = softThresh( in, thresh )
  %
  % Inputs:
  % in - scalar or array
  % thresh - the soft thresholding value
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  out = sign(in) .* max( ( abs(in) - thresh ), 0 );
  
end
