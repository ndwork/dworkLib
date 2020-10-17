
function out = mri_gridRecon( data, traj, sImg, varargin )
  % out = mri_gridRecon( data, traj, sImg [, 'alpha', alpha, 'W', W, 'nC', nC, ...
  %   'weights', weights ] )
  %
  % data - an NxC array of data values where C is the number of
  %   coils
  %
  % Outputs:
  % out - an sImg x C array that is the set of reconstructions
  %
  % Written by Nicholas Dwork, Copyright 2020
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  defaultAlpha = 1.5;
  defaultW = 8;
  defaultNc = 500;

  checknum = @(x) isnumeric(x) && isscalar(x) && (x >= 1);
  p = inputParser;
  p.addParameter( 'alpha', defaultAlpha, @(x) numel(x) == 0 || checknum(x) );
  p.addParameter( 'W', defaultW, @(x) numel(x) == 0 || checknum(x) );
  p.addParameter( 'nC', defaultNc, @(x) numel(x) == 0 || checknum(x) );
  p.addParameter( 'weights', [], @isnumeric );
  p.addParameter( 'verbose', false, @(x) islogical(x) || isnumeric(x) );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;
  weights = p.Results.weights;
  verbose = p.Results.verbose;

  if numel( alpha ) == 0, alpha = defaultAlpha; end
  if numel( W ) == 0, W = defaultW; end
  if numel( nC ) == 0, nC = defaultNc; end

  nCoils = size( data, 2 );

  if numel( weights ) == 0
    weights = makePrecompWeights_2D( traj, sImg, ...
      'alpha', alpha, 'W', W, 'nC', nC, 'verbose', verbose );
  end

  coilRecons = cell( 1, 1, nCoils );
  parfor gCoilIndx = 1 : nCoils
    coilData = data(:,gCoilIndx);
    thisRecon = grid_2D( coilData(:), traj, sImg, weights, ...
      'alpha', alpha, 'W', W, 'nC', nC );
    coilRecons{1,1,gCoilIndx} = thisRecon;
  end
  out = cell2mat( coilRecons );

end