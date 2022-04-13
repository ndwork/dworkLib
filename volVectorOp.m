

function out = volVectorOp( volume, vector, varargin )
  % This function performs a volume-vector operation
  %
  % out = volVectorOp( volume, vector [, dim, 'op', op ] )
  %
  % Inputs:
  % volume - an array
  % vector - a 1D array with number of elements equal to dimension of
  %   vector product
  %
  % Optional Input:
  % dim - Dimension along which to perform the volume-vector product
  %   (default is last dimension)
  % op - The operation to be performed
  %      Options: 'multiplication', 'subtraction', 'power'
  %
  % Written by Nicholas Dwork, Copyright 2017
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addParameter( 'dim', ndims(volume), @positive );
  p.addParameter( 'op', 'multiplication', @(x) true );
  p.parse( varargin{:} );
  dim = p.Results.dim;
  op = p.Results.op;

  newShape = ones( 1, ndims(volume) );
  newShape(dim) = numel(vector);
  newVec = reshape( vector, newShape);

  switch op
    case 'subtraction'
      out = bsxfun( @minus, volume, newVec );
    case 'multiplication'
      out = bsxfun( @times, volume, newVec );
    case 'power'
      out = bsxfun( @power, volume, newVec );
    otherwise
      error('Operation wasn''t recognized');
  end
end
