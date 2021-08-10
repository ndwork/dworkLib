
function recons = mri_gridRecon( kData, trajs, sImg, varargin )
  % recons = mri_gridRecon( kData, trajs, sImg [, 'alpha', alpha, 'W', W, 'nC', nC, ...
  %   'weights', weights, 'verbose', verbose ] )
  %
  % Performs a gridding recon for each coil
  %
  % Inputs:
  % kData - A complex N x C x D1 x D2 x ... x DN array of data values
  % trajs - a vector of N x 2 x C x D1 x D2 x ... x DN specifying the k-space trajectory
  %   If trajs is size N x 2, it is assumed the same trajectory was used for all
  %     other dimensions.
  %   C is the number of coils used during the acquisition.
  % sImg - the size of the output image
  %
  % Optional Inputs:
  % weights - the density compensation weights to be used
  %
  % Output:
  % recons is an array of reconstructed images of size sImg(1) x sImg(2) x C x D1 x D2 x ... x DN
  %
  % Written by Nicholas Dwork - Copyright 2018
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  recons = [];
  error( 'This function has been renamed to mri_reconGridding' );
end
