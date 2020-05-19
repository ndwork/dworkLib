
function out = normL2L1( in, varargin )
  % out = normL2L1( in, [ 'ws', ws ] )
  %
  % Calculates the L_{2,1} norm of the input.  It is assumed that
  % the last dimension of in represents the groups
  %
  % Inputs:
  % in: a multi-dimensional array.  The last dimension of in represents
  %     the groups.  That is, the L2 norm is calculated on the last dimension
  %
  % Optional Inputs:
  % ws: weights for a weighted L1 norm
  %
  % Outputs:
  % out: a scalar value greater than or equal to 0
  %
  % Written by Nicholas Dwork - Copyright 2020
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addParameter( 'ws', 1, @isnumeric );
  p.parse( varargin{:} );
  ws = p.Results.ws;
  
  normsL2 = norms( in, 2, ndims(in) );
  out = sum( ws(:) .* normsL2(:) );
end
