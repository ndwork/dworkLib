
function out = iGridT( k, traj, N, varargin )
  % out = iGridT( k, traj, N, [ 'alpha', alpha, 'w', w, 'nC', nC ] )
  %
  % Gridding (without density correction) is the adjoint of MRI encoding
  % (often called inverse gridding).  This function applies the transpose 
  % of inverse gridding to the input data.
  % Detailed in the following document http://nicholasdwork.com/tutorials/dworkGridding.pdf
  %
  % Inputs:
  % k is a 1D array of M elements specifying the k-space data values
  % traj is a MxV array specifying the k-space trajectory.
  %   V is the number of dimensions of the trajectory elements
  %   The first/second/third column is kx/ky/kz
  %   The units are normalized to [-0.5,0.5).
  % N is the size of the output image
  %   If N is a scalar, then the final image is assumed to be square
  %   If N is a 2 element array, then N = [Nx Ny Nz]
  % alpha is an optional float parameter specifying the oversampling factor
  % w is an integer specifying the kernel's width
  % nC specifies the number of samples in the kernel
  %
  % Written by Nicholas Dwork, Copyright 2015
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if size( traj, 1 ) == 1, traj=transpose(traj); end
  
  nD = size( traj, 2 );
  if nD == 1
    if size( k, 1 ) == 1, k=transpose(k); end
    out = iGridT_1D( k, traj, N, varargin{:} );
  elseif nD == 2
    out = iGridT_2D( k, traj, N, varargin{:} );
  elseif nD == 3
    out = iGridT_3D( k, traj, N, varargin{:} );
  end
end
