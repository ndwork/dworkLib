
function [recon,oValues] = csReconFISTA( samples, lambda, varargin )
  % recon = csReconFISTA( samples, lambda [, 'nIter', nIter, ...
  %   'printEvery', printEvery, 'transformType', transformType, 'verbose', verbose, ...
  %   'wavSplit', wavSpit ] )
  %
  % This routine minimizes 0.5 * || Ax - b ||_2^2 + lambda || W x ||_1
  %   where A is sampleMask * Fourier Transform * real part, and
  %   W is the wavelet transform.
  %
  % Inputs:
  % samples - a 2D array that is zero wherever a sample wasn't acquired
  % lambda - regularization parameter
  %
  % Optional Inputs:
  % debug - if true, reduces the default number of iterations to 30 and forces verbose
  %         statements during optimization
  % nIter - the number of iterations that FISTA will perform (default is 100)
  % printEvery - FISTA prints a verbose statement every printEvery iterations
  % transformType - either 'curvelet', 'wavelet', or 'wavCurv'.
  % verbose - if true, prints informative statements
  % waveletType - (deprecated) either 'Daubechies-4' (default) or 'Haar'
  %   If transformType is set, it overrides this option.
  %
  % Written by Nicholas Dwork - Copyright 2017
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  wavSplit = zeros(4);  wavSplit(1,1) = 1;
  %wavSplit = [1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0;];
  %wavSplit = [ 1 0; 0 0; ];

  p = inputParser;
  p.addParameter( 'alg', 'fista_wLS', @(x) true );
  p.addParameter( 'checkAdjoints', false, @islogical );
  p.addParameter( 'nIter', [], @ispositive );
  p.addParameter( 'printEvery', 1, @ispositive );
  p.addParameter( 'stepSize', 1, @ispositive );
  p.addParameter( 'transformType', 'wavCurv', @(x) true );
  p.addParameter( 'wavSplit', wavSplit, @isnumeric );
  p.addParameter( 'verbose', false, @(x) isnumeric(x) || islogical(x) );
  p.addParameter( 'waveletType', [], @(x) true );
  p.parse( varargin{:} );
  alg = p.Results.alg;
  checkAdjoints = p.Results.checkAdjoints;
  nIter = p.Results.nIter;
  printEvery = p.Results.printEvery;
  stepSize = p.Results.stepSize;
  transformType = p.Results.transformType;
  waveletType = p.Results.waveletType;
  verbose = p.Results.verbose;

  if numel( transformType ) == 0, transformType = 'wavCurv'; end
  if numel( waveletType ) == 0, waveletType = 'Daubechies-4'; end

  M = ( samples ~= 0 );
  b = samples( M == 1 );
  %x0 = zeros( size( samples ) );
  x0 = FH( samples );

  function out = F( x )
    out = fftshift( fftshift( ufft2( ifftshift( ifftshift( x, 1 ), 2 ) ), 1 ), 2 );
  end

  function out = FH( y )
    out = fftshift( fftshift( uifft2( ifftshift( ifftshift( y, 1 ), 2 ) ), 1 ), 2 );
  end

  sImg = size( samples );
  nImg = prod( sImg );

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
      cx = fdct_wrapping( x, false );
      out = curvCell2Vec( cx );
    else
      cx = vec2CurvCell( x, curvCellSizes );
      out = ifdct_wrapping( cx, false );
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

  if checkAdjoints == true
    x1 = rand( size(samples) ) + 1i * rand( size( samples ) );
    if checkAdjoint( x1, sparsifier, sparsifierH ) ~= true, error( 'sparsifier adjoint is incorrect' ); end
    if checkAdjoint( x1, @F, @FH ) ~= true, error( 'FH is not the transpose of F' ); end

    xA = sparsifier( double( imread( 'cameraman.png' ) ) / 255. );
    if checkAdjoint( xA, @A, @Aadj ) ~= true, error( 'Aadj is not the transpose of A' ); end
  end


  function out = g( x )
    diff = A( x ) - b;
    out = 0.5 * norm( diff(:), 2 ).^2;
  end

  Aadjb = Aadj( b );
  function out = gGrad( x )
    out = Aadj( A( x ) ) - Aadjb;
  end

  PsiX0 = sparsifier( x0 );  % Psi is the sparsifying transformation
  nPsiX = numel( PsiX0 );
  if numel( lambda ) == 0
    lambda = nPsiX * findFractionAboveValue( abs( PsiX0(:) ), 0.05 );
  end

  proxth = @(x,t) proxL1Complex( x, t * lambda / nPsiX );

  function out = h( x )
    out = sum( abs( x(:) ) ) * lambda / nPsiX;
  end

  t = stepSize;
  if nargout > 1
    if strcmp( alg, 'pogm' )
      [xStar,oValues] = pogm( PsiX0, @gGrad, proxth, nIter, 'g', @g, 'h', @h, 't', t, ...
        'verbose', verbose, 'printEvery', printEvery );
    elseif strcmp( alg, 'fista' )
      [xStar,oValues] = fista( PsiX0, @gGrad, proxth, 'N', nIter, 'g', @g, 'h', @h, 'verbose', verbose );
    elseif strcmp( alg, 'fista_wLS' )
      [xStar,oValues] = fista_wLS( PsiX0, @g, @gGrad, proxth, 'h', @h, ...
        't0', t, 'N', nIter, 'verbose', true, 'printEvery', printEvery );
    else
      error( 'Unrecognized algorithm' );
    end

  else
    if strcmp( alg, 'pogm' )
      xStar = pogm( PsiX0, @gGrad, proxth, nIter, 't', t );
    elseif strcmp( alg, 'fista' )
      xStar = fista( PsiX0, @gGrad, proxth, 'N', nIter );
    elseif strcmp( alg, 'fista_wLS' )
      xStar = fista_wLS( PsiX0, @g, @gGrad, proxth, 't0', t, 'N', nIter, ...
        'verbose', verbose, 'printEvery', printEvery );
    else
      error( 'Unrecognized algorithm' );
    end
  end

  recon = reshape( sparsifierH( xStar ), sImg );
end
