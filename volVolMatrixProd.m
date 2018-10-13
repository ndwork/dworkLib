
function out = volVolMatrixProd( A, B )
  % out = volVolMatrixProd( A, B )
  %
  % Matrix multiplies each slice of A by each slice of B
  %
  % Inputs:
  % A - Either a 3D array or a 2D array.
  %     If a 2D array, each row is considered a matrix multiply
  % B - Either a 3D array or a 2D array.
  %     If a 2D array, each column is considered a matrix multiply
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

  sA = size(A);
  if ismatrix(A)
    A = reshape(A,[1 sA(1) sA(2)]);
    sA = size(A);
  end

  sB = size(B);
  if ismatrix(B)
    B = reshape(B,[sB(1) 1 sB(2)]);
    sB = size(B);
  end

  if sA(3) ~= sB(3)
    error('A and B must have same size for 3rd dimension');
  end;

  out = zeros( size(A,1), size(B,2), size(A,3) );
  for j = 1:size(A,2)
    out = out + bsxfun( @times, A(:,j,:), B(j,:,:) );
  end

end

