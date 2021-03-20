
function out = gridT_2D( in, traj, weights, varargin )
  % out = gridT_2D( in, traj, weights, [ 'alpha', alpha, 'W', W, 'nC', nC ] )
  %
  % Adjoint of gridding operation.  Definitions and details according to
  % http://nicholasdwork.com/tutorials/dworkGridding.pdf
  %
  % Inputs:
  %   in is a 1D array representing the Fourier values
  %   traj is a Mx2 element array specifying the k-space trajectory.
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
  %   out is the uniformly spaced data in the space domain
  %
  % Written by Nicholas Dwork - Copyright 2015
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  tmp = iGrid_2D( in, traj, varargin{:} );

  out = bsxfun( @times, tmp, weights );
end
