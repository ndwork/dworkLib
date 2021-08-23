
function out = circConvT_2D( img, kernel )
  % out = circConvT_2D( img, kernel )
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

  [M,N] = size( kernel );
  hM = floor( M / 2 );
  hN = floor( N / 2 );

  out = zeros( size( img ) );
  for col = 1 : size(img,2)
    colShift = -(col-1) + hN;

    for row = 1 : size( img, 1 )
      rowShift = -(row-1) + hM;

      shiftedOut = circshift( out, [ rowShift, colShift ] );

      shiftedOut( 1 : M, 1 : N ) = shiftedOut( 1 : M, 1 : N ) + ...
        img( row, col ) * kernel;

      out = circshift( shiftedOut, [ -rowShift, -colShift ] );
    end
  end

end
