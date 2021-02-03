
function recons = mri_gridRecon( kData, trajs, sImg, varargin )
  % recons = mri_gridRecon( kData, trajs, sImg [, 'alpha', alpha, 'W', W, 'nC', nC, ...
  %   'verbose', verbose ] )
  %
  % Performs a gridding recon for each coil
  %
  % Inputs:
  % data - Several options:
  %   A complex N x D1 x D2 x ... x DN array of data values
  % trajs - a vector of N x 2 x D1 x D2 x ... x DN specifying the k-space trajectory
  %   If trajs is size N x 2, it is assumed the same trajectory was used for all
  %   other dimensions.
  % sImg - the size of the output image
  %
  % Output:
  % recons is an array of reconstructed images of size sImg(1) x sImg(2) x D1 x D2 x ... x DN
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
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;
  verbose = p.Results.verbose;

  sData = size( kData );
  nImgs = prod( sData(2:end) );
  kData = reshape( kData, [ sData(1), nImgs ] );
  if ismatrix( trajs )
    weights = makePrecompWeights_2D( trajs, sImg, 'alpha', alpha, 'W', W, 'nC', nC );
    trajs = repmat( trajs, [ 1 1 nImgs ] );
    weights = repmat( weights, [ 1 nImgs ] );
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

  recons = cell( 1, 1, nImgs );
  p = parforProgress( nImgs );
  parfor imgIndx = 1 : nImgs
    if verbose == true, p.progress( imgIndx ); end   %#ok<PFBNS>
    recons{imgIndx} = grid_2D( F, trajs(:,:,imgIndx), sImg, weights(:,:,imgIndx), ...
      'alpha', alpha, 'W', W, 'nC', nC );
  end
  p.clean;
  recons = cell2mat( recons );

end
