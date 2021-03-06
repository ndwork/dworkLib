
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

  if nargin < 1
    disp( 'Usage:  recons = mri_gridRecon( kData, trajs, sImg [, ''alpha'', alpha, ...' );
    disp( '          ''W'', W, ''nC'', nC, ''verbose'', verbose ' );
    return;
  end

  p = inputParser;
  p.addParameter( 'alpha', [], @ispositive );
  p.addParameter( 'W', [], @ispositive );
  p.addParameter( 'nC', [], @ispositive );
  p.addParameter( 'verbose', true, @(x) islogical(x) || isnumeric(x) );
  p.addParameter( 'weights', [], @isnumeric );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;
  verbose = p.Results.verbose;
  weights = p.Results.weights;

  sData = size( kData );
  nCoils = size( kData, 2 );
  nImgs = prod( sData(3:end) );
  kData = reshape( kData, [ sData(1), nCoils, nImgs ] );

  if ~isreal( trajs )
    trajs = [ real( trajs(:) ) imag( trajs(:) ) ];
  end
  
  if numel( weights ) == 0
    if ismatrix( trajs )
      weights = makePrecompWeights_2D( trajs, sImg, 'alpha', alpha, 'W', W, 'nC', nC );
    else
      trajs = reshape( trajs, [ sData(1) 2 nImgs ] );
      weights = cell( 1, nImgs );
      if verbose == true
        disp( 'mri_gridRecon: creating gridding weights' );
      end
      parfor imgIndx = 1 : nImgs
        weights{imgIndx} = makePrecompWeights_2D( trajs(:,:,imgIndx), sImg, ...
          'alpha', alpha, 'W', W, 'nC', nC );
      end
      weights = cell2mat( weights );
    end
  end

  if ismatrix( trajs )
    trajs = repmat( trajs, [ 1 1 nImgs ] );
  end
  if size( weights, 2 ) ~= nImgs
    weights = repmat( weights, [ 1 nImgs ] );
  end

  if nImgs > 1
    recons = cell( 1, 1, 1, nImgs );
    p = parforProgress( nImgs );
    parfor imgIndx = 1 : nImgs
      if verbose == true, p.progress( imgIndx ); end   %#ok<PFBNS>
      thisTraj = trajs(:,:,imgIndx);
      theseWeights = weights(:,imgIndx);
      recons{ 1, 1, 1, imgIndx } = grid_2D( kData(:,:,imgIndx), thisTraj, sImg, ...
          theseWeights, 'alpha', alpha, 'W', W, 'nC', nC );
    end
    recons = cell2mat( recons );
    p.clean;

  else
    recons = grid_2D( kData, trajs, sImg, weights, 'alpha', alpha, 'W', W, 'nC', nC );

  end

end
