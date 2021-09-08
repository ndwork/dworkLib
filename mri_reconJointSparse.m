
function [recons,objectiveValues] = mri_reconJointSparse( data, traj, sImg, ...
  lambda, varargin )
  % recons = mri_reconJointSparse( data, traj, sImg, lambda [, 
  %   'alpha', alpha, 'W', W, 'nC', nC ] )
  %
  % Inputs: 
  % data - an NxC array of data values where N is the number of data points and
  %   C is the number of coils
  % traj - a complex vector of N length specifying the k-space trajectory
  % sImg - the size of the output image
  % lambda - regularization parameter
  %
  % Optional Inputs:
  % alpha - oversampling factor in Gridding
  % W - window width (in pixels)
  % nC - number of points to use in convolution kernel
  %
  % Written by Nicholas Dwork, Copyright 2020
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addParameter( 'alpha', [], @ispositive );
  p.addParameter( 'W', [], @ispositive );
  p.addParameter( 'nC', [], @ispositive );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;

  nTraj = size( traj, 1 );
  nCoils = size( data, 2 );
  if size( data, 1 ) ~= nTraj, error( 'Wrong input dimensions' ); end

  nPixels = prod( sImg );
  lambdaScaled = lambda / ( nPixels * nCoils );

  nData = numel( data );
  function out = applyA( in, op )
    if nargin < 2 || strcmp( op, 'notransp' )
      y = reshape( in, [ sImg nCoils ] );
      out = cell( 1, nCoils );
      parfor cIndx = 1 : nCoils
        thisImg = y( :, :, cIndx );
        out{1,cIndx} = iGrid( thisImg, traj, 'alpha', alpha, 'W', W, 'nC', nC );
      end
      out = cell2mat( out );

    elseif nargin > 1 && strcmp( op, 'transp' )
      kData = reshape( in, [ nTraj nCoils ] );
      out = cell( 1, 1, nCoils );
      parfor cIndx = 1 : nCoils
        out{1,1,cIndx} = iGridT( kData(:,cIndx), traj, sImg, 'alpha', alpha, 'W', W, 'nC', nC );
      end
      out = cell2mat( out );

    else
      error( 'Unrecognized op' );
    end

    out = out(:) / sqrt( nData );
  end

  doCheckAdjoint = false;
  if doCheckAdjoint == true
    tmp = zeros( [ sImg nCoils ] );
    [checkA,adjErrA] = checkAdjoint( tmp, @applyA );
    if checkA ~= true
      error([ 'applyA failed adjoint test with error ', num2str(adjErrA) ]);
    end
  end

  b = data(:) / sqrt( nData );
  function out = g( in )
    Ain = applyA( in );
    out = 0.5 * norm( Ain(:) - b ).^2;
  end

  ATb = applyA( b, 'transp' );
  function out = gGrad( in )
    Ain = applyA( in );
    ATAin = applyA( Ain, 'transp' );  
    out = ATAin - ATb;
  end

  function out = h( in )
    X = reshape( in, [ sImg nCoils ] );
    WX = zeros( [ sImg nCoils ] );
    parfor cIndx = 1 : nCoils
      WX(:,:,cIndx) = wtDaubechies2( X(:,:,cIndx), wavSplit );
    end
    out = lambdaScaled * normL2L1( WX );
  end

  wavSplit = zeros(4);  wavSplit(1,1) = 1;
  function out = proxth( in, t )
    X = reshape( in, [ sImg nCoils ] );
    WX = zeros( [ sImg nCoils ] );
    parfor cIndx = 1 : nCoils
      WX(:,:,cIndx) = wtDaubechies2( X(:,:,cIndx), wavSplit );
    end
    proxL2L1WX = proxL2L1( WX, t * lambda / prod( [ sImg nCoils ] ) );
    out = zeros( [ sImg nCoils ] );
    parfor cIndx = 1 : nCoils
      out(:,:,cIndx) = iwtDaubechies2( proxL2L1WX(:,:,cIndx), wavSplit );
    end
    out = out( : );
  end

  % Setup inputs to fista_wLS
  innerProd = @(x,y) real( dotP( x, y ) );
  initializeWithGridding = true;
  if initializeWithGridding == false
    x0 = zeros( [ sImg nCoils ] );
  else
    x0 = mri_gridRecon( data, traj, sImg, 'alpha', alpha, 'W', W, 'nC', nC );
  end

  % minimize || A y - b ||_2^2 + lambda || y ||_{2,1}
  % Equivalently, minimize g(y) + h(y) where
  % g(y) = || A y - b ||_2^2  and  h(y) = lambda || y ||_{2,1}

  ATA = @(x) applyA( applyA( x ), 'transp' );
  L = powerIteration( ATA, ones( size( x0 ) ), 'symmetric', true );
  minStep = 0.95 / L;
  t0 = 10 / L;
  
  nIter = 30;
  [xStar,objectiveValues] = fista_wLS( x0(:), @g, @gGrad, @proxth, 'N', nIter, ...
    't0', t0, 'innerProd', innerProd, 'minStep', minStep, 'h', @h, 'verbose', true );

  recons = reshape( xStar, [ sImg nCoils ] );
end
