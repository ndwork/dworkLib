
function recons = mri_reconGridding( kData, traj, sImg, varargin )
  % recons = mri_reconGridding( kData, trajs, sImg [, 'alg', alg, 'alpha', alpha, ...
  %   'W', W, 'nC', nC, 'weights', weights, 'verbose', verbose ] )
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
  % alg - algorithm for makePrecompWeights_2D
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

  if nargin < 1
    disp( 'Usage:  recons = mri_reconGridding( kData, trajs, sImg [, ''alpha'', alpha, ...' );
    disp( '          ''W'', W, ''nC'', nC, ''verbose'', verbose ' );
    return;
  end

  p = inputParser;
  p.addParameter( 'alg', [], @(x) true );
  p.addParameter( 'alpha', [], @ispositive );
  p.addParameter( 'W', [], @ispositive );
  p.addParameter( 'nC', [], @ispositive );
  p.addParameter( 'verbose', true, @(x) islogical(x) || isnumeric(x) );
  p.addParameter( 'weights', [], @isnumeric );
  p.parse( varargin{:} );
  alg = p.Results.alg;
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;
  verbose = p.Results.verbose;
  weights = p.Results.weights;

  sData = size( kData );
  nCoils = sData(2);
  nImgs = prod( sData(3:end) );
  kData = reshape( kData, [ sData(1) nCoils nImgs ] );

  if ~isreal( traj )
    sTraj = size( traj );
    traj = [ real( traj(:) ) imag( traj(:) ) ];
    traj = reshape( traj, [ sTraj 2 ] );
    traj = permute( traj, [ 1 ndims(traj) 2:ndims(traj)-1 ] );
  end

  if numel( weights ) == 0
    if ismatrix( traj )
      weights = makePrecompWeights_2D( traj, 'sImg', sImg, 'alg', alg, ...
        'alpha', alpha, 'W', W, 'nC', nC );
    else
      traj = reshape( traj, [ sData(1) 2 nImgs ] );
      weights = cell( 1, nImgs );
      if verbose == true
        disp( 'mri_gridRecon: creating gridding weights' );
      end
      parfor imgIndx = 1 : nImgs
        weights{imgIndx} = makePrecompWeights_2D( traj(:,:,imgIndx), 'sImg', sImg, ...
          'alg', alg, 'alpha', alpha, 'W', W, 'nC', nC );
      end
      weights = cell2mat( weights );
    end
  end

  if ismatrix( traj )
    recons = grid_2D( kData, traj, sImg, weights, 'alpha', alpha, 'W', W, 'nC', nC );

  else
    nImgs = size( traj, 3 );
    recons = cell( 1, 1, 1, nImgs );
    p = parforProgress( nImgs );
    parfor imgIndx = 1 : nImgs
      if verbose == true, p.progress( imgIndx ); end   %#ok<PFBNS>
      thisTraj = traj(:,:,imgIndx);
      theseWeights = weights(:,imgIndx);
      recons{ 1, 1, 1, imgIndx } = grid_2D( kData(:,:,imgIndx), thisTraj, sImg, ...
          theseWeights, 'alpha', alpha, 'W', W, 'nC', nC );
    end
    recons = cell2mat( recons );
    p.clean;
  end

  recons = reshape( recons, [ sImg sData(2:end) ] );
end
