
function out = huber( in, mu )
  % out = huber( in, mu )
  %
  % Calculates the Huber penalty loss of in
  %
  % Inputs:
  % in - a 1D array
  % mu - the scalar Huber parameter
  %
  % Written by Nicholas Dwork - Copyright 2021
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  if nargin < 1
    disp( 'Usage:  out = huber( in, mu )' );
    out = [];
    return;
  end

  out = in .* in * 0.5 / mu;
  out( abs( in ) > mu ) = abs( in( abs( in ) > mu ) ) - 0.5 * mu;
  out = sum( out(:) );

end
