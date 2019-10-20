
function out = add2Diag( A, v )
  % out = add2Diag( A, v )
  %
  % Adds vector v to the diagonal elements of matrix A
  %
  % Inputs:
  % A - a matrix of size M x N
  % v - either a scalar or a vector with min(M,N) elements in it
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  sA = size( A );

  i = 1 : min( sA );
  if numel( v ) > 1 && numel( i ) ~= numel( v )
    error( 'v has the wrong number of elements' );
  end

  indxs = sub2ind( sA, i, i );

  out = A;
  out( indxs ) = out( indxs ) + v;
end
