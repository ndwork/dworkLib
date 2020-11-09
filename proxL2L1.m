
function out = proxL2L1( in, t, weights )
  % out = proxL2L1( in, t [, weights ] )
  %
  % Returns the proximal operator of f(x) = t * weights .* L2L1( x ), where
  %   L2L1 is the (possibly weighted) L2,L1 norm.
  %
  % Inputs:
  % in - three dimensional array
  % thresh - the thresholding value
  %
  % Optional Inputs:
  % weights - the weights of the norm (should it be weighted)
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin > 2, t = t .* weights; end

  nDimsIn = ndims( in );
  normsIn = norms( in, 2, nDimsIn );

  scalingFactors = t ./ normsIn;
  scalingFactors( normsIn <= t ) = 1;

  projsOntoL2Ball = bsxfun( @times, in, scalingFactors );

  out = in - projsOntoL2Ball;
end
