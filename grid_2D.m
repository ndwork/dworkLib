
function recon = grid_2D( F, kTraj, N, varargin )
  % recon = grid_2D( F, kTraj, N, [ weights, 'alpha', alpha, 'W', W, 'nC', nC ] )
  %
  % The gridding non-uniform FFT algorithm based on EE369C notes by John Pauly
  % and Beatty et. al., IEEE TMI, 2005.  Definitions and details according to
  % http://nicholasdwork.com/tutorials/dworkGridding.pdf
  %
  % Inputs:
  %   F is a 1D array representing the Fourier values
  %   kTraj is a Mx2 array specifying the k-space trajectory or an M element complex array.
  %     If Mx2, then the first/second column is kx/ky.  Otherwise, kx/ky is real/imag.
  %     The units are normalized to [-0.5,0.5).
  %   N is a 2 element array [Ny Nx] representing the number of grid points
  %
  % Optional Inputs:
  %   weights is a 1D array; it is the pre-density compensation weights and
  %     can be generated using makePrecompWeights_2D.  Alternatively, they
  %     can be determined analytically for some sequences.
  %     By default, weights are determined with makePrecompWeights_2D.
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

  if nargin < 1
    disp( 'Usage:  recon = grid_2D( F, kTraj, N, [ weights, ''alpha'', alpha, ''W'', W, ''nC'', nC ] )' );
    if nargout > 0, recon = []; end
    return
  end

  p = inputParser;
  p.addOptional( 'weights', [], @isnumeric );
  p.addParameter( 'alpha', [], @(x) numel(x) == 0  ||  x >= 1 );
  p.addParameter( 'nC', [], @(x) numel(x) == 0  ||  ispositive(x) );
  p.addParameter( 'W', [], @(x) numel(x) == 0  ||  ispositive(x) );
  p.parse( varargin{:} );
  weights = p.Results.weights;
  alpha = p.Results.alpha;
  nC = p.Results.nC;
  W = p.Results.W;

  if numel( weights ) == 0
    weights = makePrecompWeights_2D( kTraj, 'alpha', alpha, 'W', W, 'nC', nC, 'sImg', N );
  end

  if numel( weights ) > 1 || weight ~= 1
    F = bsxfun( @times, F, weights );
  end

  recon = iGridT_2D( F, kTraj, N, 'alpha', alpha, 'W', W, 'nC', nC );
end
