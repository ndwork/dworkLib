
function [recon,oValues,lambda] = csReconLASSO( samples, varargin )
  % recon = csReconLASSO( samples [, 'lambda', lambda, 'nIter', nIter, 'nReweightIter', nReweightIter, ...
  %   'printEvery', printEvery, 'transformType', transformType, 'verbose', verbose, ...
  %   'waveletType', waveletType, 'w', w, 'wavSplit', wavSpit ] )
  %
  % This routine uses proximal gradient methods to minimize
  %   0.5 * || A y - b ||_2^2 + lambda || y ||_{w,1}
  %   where A is sampleMask * Fourier Transform * adjoint( Psi ).
  %   Here, Psi is either a wavelet or curvelet transform.
  %   The reconstruction returned is adjoint( Psi ) * y.
  %
  % Note that || y ||_{w,1} = w1 |y1| + w2 |y2| + ... + wN |yN|
  %
  % Iterative reweighting is done according to "Enhancing sparsity by reweighted L1
  %   minimization" by Candes, Watkin, and Boyd with the nReweightIter optional parameter.
  %
  % Inputs:
  % samples - a 2D array that is zero wherever a sample wasn't acquired
  % lambda - regularization parameter
  %
  % Optional Inputs:
  % alg - choose the algorithm to perform the minimization.  Default is 'FISTA_wLS'.
  %   Options are FISTA_wLS.
  % nIter - the number of iterations that FISTA will perform (default is 100)
  % nReweightIter - number of reweighting iterations
  % printEvery - FISTA prints a verbose statement every printEvery iterations
  % transformType - Choose the sparifying transformation.  Default is WavCurv (which uses
  %   the wavelet transform and the curvelet transform as a redundant dictionary).
  %   Options are:  'curvlet', 'wavelet', or 'wavCurv'.
  %     curvelet - Discrete curvelet transform (requires CurveLab; see curvelet.org)
  %     wavelet - Wavelet transform (type specified by waveletType parameter)
  %     wavCurv - Redundant dictionary of wavelet and curvelet (default)
  % verbose - if true, prints informative statements
  % waveletType - either 'Daubechies-4' (default) or 'Haar'
  %
  % Written by Nicholas Dwork - Copyright 2017
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  wavSplit = zeros(4);  wavSplit(1,1) = 1;

  p = inputParser;
  p.addParameter( 'lambda', [], @isnumeric );
  p.addParameter( 'alg', 'fista_wLS', @(x) true );
  p.addParameter( 'checkAdjoints', false, @islogical );
  p.addParameter( 'nIter', [], @ispositive );
  p.addParameter( 'nReweightIter', [], @ispositive );
  p.addParameter( 'printEvery', 1, @ispositive );
  p.addParameter( 'stepSize', [], @ispositive );
  p.addParameter( 'transformType', 'wavCurv', @(x) true );
  p.addParameter( 'verbose', false, @(x) isnumeric(x) || islogical(x) );
  p.addParameter( 'w', 1, @isnumeric );
  p.addParameter( 'waveletType', [], @(x) true );
  p.addParameter( 'wavSplit', wavSplit, @isnumeric );
  p.parse( varargin{:} );
  lambda = p.Results.lambda;
  alg = p.Results.alg;
  checkAdjoints = p.Results.checkAdjoints;
  nIter = p.Results.nIter;
  nReweightIter = p.Results.nReweightIter;
  printEvery = p.Results.printEvery;
  stepSize = p.Results.stepSize;
  transformType = p.Results.transformType;
  verbose = p.Results.verbose;
  w = p.Results.w;
  waveletType = p.Results.waveletType;
  wavSplit = p.Results.wavSplit;

  if numel( nReweightIter ) == 0 || nReweightIter == 0
    if numel( lambda ) == 0 || lambda == 0
      nReweightIter = 4;
    else
      nReweightIter = 1;
    end
  end
  reweightEpsilon = 0.1;

  if numel( transformType ) == 0, transformType = 'wavCurv'; end
  if numel( waveletType ) == 0, waveletType = 'Daubechies-4'; end

  sImg = size( samples );
  nImg = prod( sImg );

  function out = F( x )
    out = fftshift2( fft2( ifftshift2( x ) ) );
  end

  function out = FH( y )
    out = fftshift2( fft2h( ifftshift2( y ) ) );
  end

  function out = findCellSizes( in )
    if iscell( in )
      out = cell( size(in) );
      for indx = 1 : numel( in )
        theseCellSizes = findCellSizes( in{ indx } );
        out{ indx } = theseCellSizes;
      end
    else
      out = size( in );
    end
  end

  function nTotal = sumCurvCellSizes( cellSizes )
    nTotal = 0;
    if iscell( cellSizes )
      for indx = 1 : numel( cellSizes )
        thisCell = cellSizes{ indx };  
        nTotal = nTotal + sumCurvCellSizes( thisCell );
      end
    else
      nTotal = nTotal + prod( cellSizes );
    end
  end

  function out = curvCell2Vec( cx )
    if iscell( cx )
      out = cell( numel( cx ), 1 );
      for indx = 1 : numel( cx )
        out{ indx } = curvCell2Vec( cx{ indx } );
      end
      out = cell2mat( out );
    else
      out = cx(:);
    end
  end

  function out = vec2CurvCell( v, curvCellSizes )
    if iscell( curvCellSizes )
      out = cell( size( curvCellSizes ) );
      thisIndx = 1;
      for indx = 1 : numel( curvCellSizes )
        nSubVec = sumCurvCellSizes( curvCellSizes{ indx } );
        subVec = v( thisIndx : thisIndx + nSubVec - 1 );
        out{ indx } = vec2CurvCell( subVec, curvCellSizes{ indx } );
        thisIndx = thisIndx + nSubVec;
      end
    else
      out = reshape( v, curvCellSizes );
    end
  end

  function out = curvelet( x, type )
    if nargin < 2 || strcmp( type, 'notransp' )
      curvCells = fdct_wrapping( x, false );
      out = curvCell2Vec( curvCells );
    else
      curvCells = vec2CurvCell( x, curvCellSizes );
      out = ifdct_wrapping( curvCells, false );
    end
  end

  if strcmp( waveletType, 'Daubechies-4' )
    wavTrans = @(x) wtDaubechies2( x, wavSplit );
    wavTransH = @(y) iwtDaubechies2( y, wavSplit );
  elseif strcmp( waveletType, 'Haar' )
    wavTrans = @(x) wtHaar2( x, wavSplit );
    wavTransH = @(y) iwtHaar2( y, wavSplit );
  else
    error( 'Unrecognized wavelet type' );
  end

  function out = wavCurv( x, type )
    if nargin < 2 || strcmp( type, 'notransp' )
      x = reshape( x, sImg );
      wx = wavTrans( x );
      cx = curvelet( x );
      out = [ wx(:); cx(:); ];
    elseif strcmp( type, 'transp' )
      x1 = reshape( x( 1 : nImg ), sImg );
      x2 = x( nImg + 1 : end );
      out = wavTransH( x1 ) + curvelet( x2, 'transp' );
      out = out(:);
    end
  end

  if strcmp( transformType, 'curvelet' ) || strcmp( transformType, 'wavCurv' )
    tmp = fdct_wrapping( FH( samples ), false );
    curvCellSizes = findCellSizes( tmp );

    if strcmp( transformType, 'wavCurv' )
      sparsifier = @(x) wavCurv( x );
      sparsifierH = @(x) wavCurv( x, 'transp' );
    else
      sparsifier = @(x) curvelet( x );
      sparsifierH = @(x) curvelet( x, 'transp' );
    end

  elseif strcmp( transformType, 'wavelet' )
    sparsifier = wavTrans;
    sparsifierH = wavTransH;

  else
    error( 'Unrecognized transform type' );
  end

  M = ( samples ~= 0 );

  function out = A( x )
    PsiHx = sparsifierH( x );
    PsiHx = reshape( PsiHx, sImg );
    FPsiHx = F( PsiHx );
    out = FPsiHx( M == 1 );
  end

  MTy = zeros( size(M) );
  function out = Aadj( y )
    MTy( M == 1 ) = y;
    FHMTy = FH( MTy );
    out = sparsifier( FHMTy );
  end

  applyATA = @(in) Aadj( A( in ) );

  if checkAdjoints == true
    innerProd = @(x,y) real( y(:)' * x(:) );

    x1 = rand( size(samples) ) + 1i * rand( size( samples ) );
    if checkAdjoint( x1, @F, @FH, 'innerProd', innerProd ) ~= true
      error( 'FH is not the transpose of F' );
    end

    if checkAdjoint( x1, sparsifier, sparsifierH, 'innerProd', innerProd ) ~= true
      error( 'sparsifierH is not the transpose of sparsifier' );
    end

    xA = sparsifier( double( imread( 'cameraman.png' ) ) / 255. );
    if checkAdjoint( xA, @A, @Aadj, 'innerProd', innerProd ) ~= true
      error( 'Aadj is not the transpose of A' );
    end
  end

  b = samples( M == 1 );
  function out = g( x )
    diff = A( x ) - b;
    out = 0.5 * norm( diff(:) ).^2;
  end

  Aadjb = Aadj( b );
  function out = gGrad( x )
    out = Aadj( A( x ) ) - Aadjb;
  end

  %img0 = zeros( size( samples ) );
  img0 = fftshift2( ifft2( ifftshift2( samples ) ) );
  PsiImg0 = sparsifier( img0 );  % Psi is the sparsifying transformation
  PsiImg0 = PsiImg0 / norm( PsiImg0(:) ) * norm( img0(:) );
  nPsi = numel( PsiImg0 );
  if numel( lambda ) == 0  ||  lambda == 0
    lambda = nPsi ./ ( abs( PsiImg0 ) + reweightEpsilon );
  end

  proxth = @(x,t) proxL1Complex( x, ( t / nPsi ) * ( w .* lambda ) );

  function out = h( x )
    out = sum( abs( x(:) .* lambda(:) .* w(:) ) ) / nPsi;
  end

  if numel( stepSize ) == 0
    [normATA,piFlag] = powerIteration( applyATA, PsiImg0, 'symmetric', true );  %#ok<ASGLU> 
    stepSize = 0.95 / normATA;
  end

  for reweightIter = 1 : nReweightIter

    if verbose == true
      disp([ 'Working on reweighting iteration ', indx2str( reweightIter, nReweightIter ), ...
        ' of ', num2str( nReweightIter ) ]);
    end

    if reweightIter > 1
      PsiImg0 = xStar;
      lambda = nPsi ./ ( abs( PsiImg0(:) ) + reweightEpsilon );
    end

    t = stepSize;
    if nargout > 1
      if strcmp( alg, 'pogm' )
        [xStar,oValues] = pogm( PsiImg0, @gGrad, proxth, nIter, 'g', @g, 'h', @h, 't', t, ...
          'tol', 1d-8, 'verbose', verbose, 'printEvery', printEvery );
      elseif strcmp( alg, 'fista' )
        [xStar,oValues] = fista( PsiImg0, @gGrad, proxth, 'N', nIter, 'g', @g, 'h', @h, 't', t, ...
          'tol', 1d-8, 'verbose', verbose, 'printEvery', printEvery );
      elseif strcmp( alg, 'fista_wLS' )
        t = t * 10;
        [xStar,oValues] = fista_wLS( PsiImg0, @g, @gGrad, proxth, 'h', @h, ...
          't0', t, 'tol', 1d-8, 'N', nIter, 'verbose', verbose, 'printEvery', printEvery );
      else
        error( 'Unrecognized algorithm' );
      end
    else
      if strcmp( alg, 'pogm' )
        xStar = pogm( PsiImg0, @gGrad, proxth, nIter, 't', t', 'tol', 1d-8 );
      elseif strcmp( alg, 'fista' )
        xStar = fista( PsiImg0, @gGrad, proxth, 'N', nIter, 't', t, 'tol', 1d-8 );
      elseif strcmp( alg, 'fista_wLS' )
        xStar = fista_wLS( PsiImg0, @g, @gGrad, proxth, 't0', t, 'tol', 1d-8, 'N', nIter, ...
          'verbose', verbose, 'printEvery', printEvery );
      else
        error( 'Unrecognized algorithm' );
      end
    end

  end

  recon = reshape( sparsifierH( xStar ), sImg );
end
