
function out = proxL2L1( in, thresh )
  % out = proxL2L1( in, thresh )
  %
  % Returns the proximal operator of f(x) = thresh * L2L1( x )
  %
  % Inputs:
  % in - three dimensional array
  % thresh - the thresholding value
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  normsIn = norms( in, 2, 3 );

  scalingFactors = thresh ./ normsIn;
  scalingFactors( abs( normsIn ) <= thresh ) = 1;
  scalingFactors = repmat( scalingFactors, [1 size(in,3)] );
  projsOntoL2Ball = in .* scalingFactors;

  out = in - projsOntoL2Ball;
end
