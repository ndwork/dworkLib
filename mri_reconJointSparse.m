
function [recons,objectiveValues] = mri_reconJointSparse( data, traj, sImg, ...
  lambda, varargin )
  % recons = mri_reconJointSparse( data, traj, sImg, lambda [, 
  %   'alpha', alpha, 'W', W, 'nC', nC ] )
  %
  % Inputs: 
  % data - an NxC array of data values where N is the number of data points and
  %   C is the number of coils
  % traj - a complex vector of N length specifying the k-space trajectory
  % sMaps - sensitivity maps
  % lambda - regularization parameter
  %
  % Optional Inputs:
  % alpha - oversampling factor in Gridding
  % W - window width (in pixels)
  % nC - number of points to use in convolution kernel
  %
  % Written by Nicholas Dwork, Copyright 2020

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

  nData = numel( data );
  wavSplit = zeros(4);  wavSplit(1,1) = 1;
  function out = applyA( in, op )
    if nargin < 2 || strcmp( op, 'notransp' )
      y = reshape( in, [ sImg nCoils ] );
      out = cell( 1, nCoils );
      parfor coilIndx = 1 : nCoils
        thisImg = y( :, :, coilIndx );
        WhImg = iwtDeaubechies2( thisImg, wavSplit );
        out{coilIndx} = iGrid( WhImg, traj, 'alpha', alpha, 'W', W, 'nC', nC );
      end
      out = cell2mat( out );

    elseif nargin > 1 && strcmp( op, 'transp' )
      kData = reshape( in, [ nTraj nCoils ] );
      WGihk = cell( 1, 1, nCoils );
      parfor coilIndx = 1 : nCoils
        Gihk = iGridT( kData(:,coilIndx), traj, sImg, 'alpha', alpha, 'W', W, 'nC', nC );
        WGihk{coilIndx} = wtDeaubechies2( Gihk, wavSplit );
      end
      out = cell2mat( WGihk );

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

  function out = proxth( in, t )
    X = reshape( in, [ sImg nCoils ] );
    out = proxL2L1( X, t * lambda );
    out = out( : );
  end

  function out = h( in )
    X = reshape( in, [ prod( sImg ) nCoils ] );
    out = normL2L1( X );
  end

  % Setup inputs to fista_wLS
  innerProd = @(x,y) real( dotP( x, y ) );
  y0 = zeros( [ sImg nCoils ] );

  % minimize || A y - b ||_2^2 + lambda || y ||_{2,1}
  % Equivalently, minimize g(y) + h(y) where
  % g(y) = || A y - b ||_2^2  and  h(y) = lambda || y ||_{2,1}

  nIter = 30;
nIter = 2;
  [yStar,objectiveValues] = fista_wLS( y0(:), @g, @gGrad, @proxth, 'N', nIter, ...
    't0', 1d-2, 'innerProd', innerProd, 'minStep', 1d-11, 'h', @h, 'verbose', true );

  yStar = reshape( yStar, [ sImg nCoils ] );
  recons = cell( 1, 1, nCoils );
  parfor coilIndx = 1 : nCoils
    recons{coilIndx} = iwtDeaubechies2( yStar(:,:,coilIndx), wavSplit );
  end
  recons = cell2mat( recons );
end
