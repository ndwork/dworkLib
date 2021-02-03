
function [recon,objValues] = mri_reconLowRankPlusJointSparse( data, trajs, sImg, ...
  lambda, sigma, varargin )
  % [recon,objValues] = mri_reconLowRankPlusJointSparse( data, trajs, sImg, ...
  %   lambda, sigma [, 'alpha', alpha, 'W', W, 'nC', nC ] )
  %
  % Inputs: 
  % data - an NxCxT array of data values where N is the number of data points and
  %   C is the number of coils
  % trajs - a vector of Nx2xT length specifying the k-space trajectory
  %   T is the number of time points
  %   If trajs is of size Nx2, it is assumed that the same trajectory was used for
  %     all time points.
  % sImg - the size of the output image
  % lambda - joint sparsity regularization parameter
  % sigma - low-rank regularization parameter
  %
  % Optional Inputs:
  % alpha - oversampling factor in Gridding
  % W - window width (in pixels) for Gridding
  % nC - number of points to use in convolution kernel in Gridding
  %
  % Outputs:
  % recon - a complex array of size sImg(1)xsImg(2)xCxT
  %
  % Written by Nicholas Dwork, Copyright 2020
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp( 'Usage:  [recon,objValues] = mri_reconLowRankPlusJointSparse( data, trajs, sImg, ...' );
    disp( '          lambda, sigma [, ''alpha'', alpha, ''W'', W, ''nC'', nC ] )' );
    return
  end

  p = inputParser;
  p.addParameter( 'nIterPDHG', 100, @ispositive );
  p.addParameter( 'alpha', [], @ispositive );
  p.addParameter( 'W', [], @ispositive );
  p.addParameter( 'nC', [], @ispositive );
  p.addParameter( 'verbose', true, @(x) islogical(x) || isnumeric(x) );
  p.parse( varargin{:} );
  nIterPDHG = p.Results.nIterPDHG;
  alpha = p.Results.alpha;
  W = p.Results.W;
  nC = p.Results.nC;
  verbose = p.Results.verbose;

  nData = numel( data );
  sData = size( data );
  nTraj = size( data, 1 );   %#ok<NASGU>
  nCoils = size( data, 2 );
  nTimes = size( data, 3 );
  sOut = [ sImg nCoils nTimes ];
  nPixels = prod( sImg );
  nOut = prod( sOut );

  if ismatrix( trajs )
    weights = sqrt( makePrecompWeights_2D( trajs, sImg, 'alpha', alpha, 'W', W, 'nC', nC ) );
    trajs = repmat( trajs, [1 1 nTimes] );
    weights = repmat( weights, [1 nTimes] );
  else
    weights = cell( 1, nTimes );
    parfor ti = 1 : nTimes
      weights{1,ti} = sqrt( makePrecompWeights_2D( trajs(:,:,ti), sImg, ...
        'alpha', alpha, 'W', W, 'nC', nC ) );
    end
    weights = cell2mat( weights );
  end

  sigmaScaled = sigma / ( nPixels * nCoils );

  function out = applyGi( in, type )
    if nargin < 2 || strcmp( type, 'notransp' )
      in = reshape( in, sOut );
      out = cell( 1, nCoils, nTimes );
      parfor timeIndx = 1 : nTimes
        traj = trajs(:,:,timeIndx);
        for coilIndx = 1 : nCoils
          iGridIn = iGrid( in(:,:,coilIndx,timeIndx), traj, 'alpha', alpha, 'W', W, 'nC', nC );
          out{1,coilIndx,timeIndx} = weights(:,timeIndx) .* iGridIn;
        end
      end
      out = cell2mat( out );

    else
      out = cell( 1, 1, nCoils, nTimes );
      in = reshape( in, sData );
      parfor timeIndx = 1 : nTimes
        traj = trajs(:,:,timeIndx);
        for coilIndx = 1 : nCoils
          wIn = in(:,coilIndx,timeIndx) .* conj( weights(:,timeIndx) );
          out{1,1,coilIndx,timeIndx} = iGridT( wIn, ...
            traj, sImg, 'alpha', alpha, 'W', W, 'nC', nC );
        end
      end
      out = cell2mat( out );
    end

    out = out / sqrt( nData );
  end

  wavSplit = zeros( 4 );  wavSplit(1,1) = 1;
  function out = applyWav( in, type )
    out = cell( 1, 1, nCoils, nTimes );
    if nargin < 2 || strcmp( type, 'notransp' )
      parfor timeIndx = 1 : nTimes
        for coilIndx = 1 : nCoils
          out{1,1,coilIndx,timeIndx} = wtDaubechies2( in(:,:,coilIndx,timeIndx), wavSplit );
        end
      end
    else
      parfor timeIndx = 1 : nTimes
        for coilIndx = 1 : nCoils
          out{1,1,coilIndx,timeIndx} = iwtDaubechies2( in(:,:,coilIndx,timeIndx), wavSplit );
        end
      end
    end
    out = cell2mat( out );
  end

  function out = applyA( in, t )
    if nargin < 2 || strcmp( t, 'notransp' )
      in = reshape( in, sOut );
      GiIn = applyGi( in );
      out = [ GiIn(:); in(:); ];
    else
      in1 = reshape( in( 1 : nData ), sData );
      in2 = reshape( in( nData+1 : nData + nOut ), sOut );
      GiTIn1 = applyGi( in1, 'transp' );
      out = GiTIn1 + in2;
    end
  end

  doCheckAdjoint = false;
  if doCheckAdjoint == true
    tmp = zeros( [ sImg nCoils nTimes ] );
    [checkA,adjErrA] = checkAdjoint( tmp, @applyA );
    if checkA ~= true
      error([ 'applyA failed adjoint test with error ', num2str(adjErrA) ]);
    end
  end

  function out = f( in )
    out = 0;
    Win = applyWav( in );
    for timeIndx = 1 : nTimes
      out = normL2L1( Win(:,:,:,timeIndx) );
    end
    out = sigmaScaled * out;
  end

  function out = proxf( in, t )
    Win = applyWav( in );
    proxWin = cell(1,1,1,nTimes);
    parfor timeIndx = 1 : nTimes
      tmp = proxL2L1( Win(:,:,:,timeIndx), t * sigmaScaled );
      proxWin{1,1,1,timeIndx} = tmp;
    end
    proxWin = cell2mat( proxWin );
    out = in - applyWav( Win - proxWin, 'transp' );
  end

  b = zeros( size( data ) );
  for cIndx = 1 : nCoils
    tmp = squeeze( data(:,cIndx,:) );
    b(:,cIndx,:) = bsxfun( @rdivide, tmp, weights );
  end
  b = b / sqrt( nData );
  b = b(:);
  function out = g( in )
    in1 = in( 1 : nData );
    in2 = reshape( in( nData + 1 : nData + nOut ), [ nPixels, nCoils, nTimes ] );

    out1 = 0.5 * norm( in1 - b )^2;

    out2 = 0;
    in2 = permute( in2, [ 1 3 2 ] );
    for coilIndx = 1 : nCoils
      xMatrix = in2(:,:,coilIndx);
      out2 = out2 + nucNorm( xMatrix );
    end
    out2 = lambda * out2;

    out = out1 + out2;
  end

  function out = proxConjG2( in, t )
    out = cell( 1, 1, nCoils );
    in = permute( in, [ 1 3 2 ] );
    parfor coilIndx = 1 : nCoils
      out{1,1,coilIndx} = proxConjNucNorm( in(:,:,coilIndx), t, lambda );
    end
    out = cell2mat( out );
    out = permute( out, [ 1 3 2 ] );
  end

  function out = proxgConj( in, t )
    in1 = in( 1 : nData );
    in2 = in ( nData + 1 : nData + nOut );

    out1 = proxConjL2Sq( in1, t, 1, b );

    in2 = reshape( in2, [ nPixels nCoils nTimes ] );
    out2 = proxConjG2( in2, t );

    out = [ out1(:); out2(:); ];
  end


  x0 = cell( 1, 1, 1, nTimes );
  parfor ti = 1 : nTimes
    x0{1,1,1,ti} = mri_gridRecon( data(:,:,ti), trajs(:,:,ti), sImg, ...
      'alpha', alpha, 'W', W, 'nC', nC );
  end
  x0 = cell2mat( x0 );
  %x0 = zeros( sOut );

  %normA = powerIteration( @applyA, ones( size( x0 ) ), 'maxIters', 100 );
  normA = 1.5;
  pdhgTau = 0.99 / ( normA * normA );  % value by Ong et al.
  theta = 1;

  adaptive = false;
  if adaptive == true
    pdhgSigma = 1 / normA;
    if nargout > 1
      [xStar,objValues] = pdhgAdaptive( x0, @proxf, @proxgConj, pdhgTau, 'sigma', pdhgSigma, ...
        'A', @applyA, 'f', @f, 'g', @g, 'N', nIterPDHG, 'normA', normA, 'verbose', verbose );
    else
      xStar = pdhgAdaptive( x0, @proxf, @proxgConj, pdhgTau, 'sigma', pdhgSigma, ...
        'A', @applyA, 'N', nIterPDHG, 'normA', normA, 'verbose', verbose );
    end
  else
    if nargout > 1
      [xStar,objValues] = pdhg( x0, @proxf, @proxgConj, pdhgTau, 'theta', theta, ...
        'A', @applyA, 'f', @f, 'g', @g, 'N', nIterPDHG, 'normA', normA, 'verbose', verbose );
    else
      [xStar,objValues] = pdhg( x0, @proxf, @proxgConj, pdhgTau, 'theta', theta, ...
        'A', @applyA, 'N', nIterPDHG, 'normA', normA, 'verbose', verbose );
    end
  end

  recon = reshape( xStar, [ sImg nCoils nTimes ] );
end


