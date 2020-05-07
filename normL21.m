
function out = normL21( in )
  % out = normL21( in )
  %
  % Calculates the L_{2,1} norm of the input.  It is assumed that
  % the last dimension of in represents the groups
  %
  % Inputs:
  % in: a multi-dimensional array.  The last dimension of in represents
  %     the groups.  That is, the L2 norm is calculated on the last dimension
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

  normsL2 = norms( in, 2, ndims(in) );
  out = sum( normsL2(:) );
end
