
function k = iGrid( img, traj, varargin )
  % k = iGrid( img, traj, [ 'alpha', alpha, 'W', W, 'nC', nC ] )
  %
  % Inverse Gridding based on Beatty et. al., IEEE TMI, 2005
  % Detailed in the following document http://nicholasdwork.com/tutorials/dworkGridding.pdf
  %
  % Inputs:
  % img is an array specifying the volume to be encoded
  % traj is a MxV array specifying the k-space trajectory.
  %   V is the number of dimensions of the img
  %   The first/second/third column is kx/ky/kz
  %   The units are normalized to [-0.5,0.5).
  % alpha is the oversampling factor [1,inf]
  % W is the window width in pixels
  % nC is the number of points to sample the convolution kernel
  %
  % Written by Nicholas Dwork - Copyright 2015
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if size( traj, 1 ) == 1, traj=transpose( traj ); end

  nD = size( traj, 2 );
  if nD == 1
    k = iGrid_1D( img, traj, varargin{:} );
  elseif nD == 2
    k = iGrid_2D( img, traj, varargin{:} );
  elseif nD == 3
    k = iGrid_3D( img, traj, varargin{:} );
  end
end
