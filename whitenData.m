
function out = whitenData( data, covMatrix )
  % out = whitenData( data, covMatrix )
  %
  % Assumes data is zero mean
  %
  % Inputs:
  % data - a multidimensional array; it is assumed that the last dimension
  %   is the random vector dimension
  % covMatrix - a covariance matrix
  %
  % Outputs:
  % out - an array of whitened vectors of size data
  %
  % Written by Nicholas Dwork - Copyright 2020
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.


  if size( data, ndims( data ) ) ~= size( covMatrix, 1 ) || ...
     size( data, ndims( data ) ) ~= size( covMatrix, 2 )
     error( 'Dimension mismatch' );
  end

  L = chol( covMatrix, 'lower' );

  sData = size( data );
  reshaped = transpose( reshape( data, [ prod( sData(1:end-1) ) sData(end) ] ) );
  tmp = L \ reshaped;

  out = reshape( tmp, sData );
end
