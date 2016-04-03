
function out = volVectorProd( volume, vector, dim )
  % This function performs a volume-vector product
  %
  % out = volVectorProd( volume, vector [, dim ] )
  %
  % Inputs:
  % volume - a 3D array of size MxNxK
  % vector - a 1D array with number of elements equal to dimension of
  %   vector product
  % Optional Input:
  % dim - Dimension along which to perform the volume-vector product
  %   (default is 3)
  %
  % Written by Nicholas Dwork

  if nargin < 3
    dim = 3;
  end
  
  if dim == 1
    newVec = reshape( vector, [numel(vector) 1 1] );
  elseif dim == 2
    newVec = reshape( vector, [1 numel(vector) 1] );
  elseif dim == 3
    newVec = reshape( vector, [1 1 numel(vector)] );
  end
  
  out = bsxfun( @times, volume, newVec );
end