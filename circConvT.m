
function out = circConvT( img, kernel )
  % out = circConvT( img, kernel )
  %
  % Calculates the adjoint (transpose) of circular convolution with kernel
  %
  % Inputs:
  % img - a 2D array
  % kernel - a small 2D array
  %
  % Outputs:
  % out - a 2D array that is the result of the adjoint of convolution with the kernel
  %
  % Written by Nicholas Dwork - Copyright 2021
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  kFlipped = flipAboutIndx( kernel, floor( size( kernel ) / 2 ) + 1 );

  out = circConv( img, kFlipped );
end
