
function out = proxL2L1( in, thresh, weights )
  % out = proxL2L1( in, thresh, weights )
  %
  % Returns the proximal operator of f(x) = thresh * L2L1( x ), where
  %   L2L1 is the (possibly weighted) L2,L1 norm.
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

  if nargin > 2
    thresh = thresh .* weights;
  end

  sIn = size( in );
  nDimsIn = numel( sIn );
  normsIn = norms( in, 2, nDimsIn );

  scalingFactors = thresh ./ normsIn;
  scalingFactors( abs( normsIn ) <= thresh ) = 1;
  repDims = ones( 1, nDimsIn );
  repDims( end ) = sIn( end );
  scalingFactors = repmat( scalingFactors, repDims );
  projsOntoL2Ball = in .* scalingFactors;

  out = in - projsOntoL2Ball;
end
