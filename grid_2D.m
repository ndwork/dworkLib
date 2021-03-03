
function recon = grid_2D( F, kTraj, N, weights, varargin )
  % recon = grid_2D( F, kTraj, N, weights, [ 'alpha', alpha, 'W', W, 'nC', nC ] )
  %
  % The gridding non-uniform fft algorithm
  % Based on EE369C notes by John Pauly and Beatty et. al., IEEE TMI, 2005
  % Adjoint of gridding operation.  Definitions and details according to
  % http://nicholasdwork.com/tutorials/dworkGridding.pdf
  %
  % Inputs:
  %   F is a 1D array representing the Fourier values
  %   kTraj is a Mx2 element array specifying the k-space trajectory.
  %     The first/second column are the kx/ky locations.
  %     The units are normalized to [-0.5,0.5).
  %   N is a 2 element array [Ny Nx] representing the number of grid points
  %   weights is a 1D array; it is the pre-density compensation weights and
  %     can be generated using makePrecompWeights_2D.  Alternatively, they
  %     can be determined analytically for some sequences.
  %
  % Optional Inputs:
  %   alpha is the oversampling factor > 1
  %   W is the window width in pixels
  %   nC is the number of points to sample the convolution kernel
  %
  % Output:
  %   recon is the uniformly spaced data in the space domain
  %
  % Written by Nicholas Dwork - Copyright 2015
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  checknum = @(x) numel(x) == 0 || ( isnumeric(x) && isscalar(x) && (x > 1) );
  p = inputParser;
  p.addParameter( 'alpha', 1.5, checknum );
  p.addParameter( 'W', [], checknum );
  p.addParameter( 'nC', [], checknum );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;

  nGrid = ceil( alpha * N );
  trueAlpha = max( nGrid ./ N );

  weightedF = bsxfun( @times, F, weights );

  padded = iGridT_2D( weightedF, kTraj, nGrid, 'alpha', trueAlpha, 'W', W, 'nC', nC );

  nOut = N(1) * N(2);
  padded = padded / ( nOut * nOut );

  recon = cropData( padded, [ N size(F,2) ] );
end

