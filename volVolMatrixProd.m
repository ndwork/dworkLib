
function out = volVolMatrixProd( A, B )
  % out = volVolMatrixProd( A, B )
  %
  % Inputs:
  % A,B - 2 3D arraysof the same sizes
  %
  % Outputs:
  % out - a 3D array of same size as A and B such that
  %  out(:,:,j) = A(:,:,j) * B(:,:,j) for all j
  %  Note that * means matrix multiplication
  %
  % Created according to http://www.alecjacobson.com/weblog/?p=4186
  %
  % Written by Nicholas Dwork - Copyright 2018
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  out = zeros( size(A) );
  for j = 1:size(A,2)
    out = out + bsxfun( @times, A(:,j,:), B(j,:,:) );
  end

end

