
function out = volVectorProd( volume, vector, dim )
  % This function performs a volume-vector product
  %
  % out = volVectorProd( volume, vector [, dim ] )
  %
  % Inputs:
  % volume - an array
  % vector - a 1D array with number of elements equal to dimension of
  %   vector product
  %
  % Optional Input:
  % dim - Dimension along which to perform the volume-vector product
  %   (default is last dimension)
  %
  % Written by Nicholas Dwork, Copyright 2017
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 3, dim=ndims(volume); end

  newShape = ones( 1, ndims(volume) );
  newShape(dim) = numel(vector);
  newVec = reshape( vector, newShape);

  out = bsxfun( @times, volume, newVec );
end