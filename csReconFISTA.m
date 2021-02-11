
function [recon,oValues] = csReconFISTA( samples, lambda, varargin )
  % recon = csReconFISTA( samples, lambda [, 'debug', debug, 'nIter', nIter, ...
  %   'printEvery', printEvery, 'verbose', verbose, 'transformType', transformType ] )
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
  % transformType - either 'Daubechies' for for Daubechies-4 (default), 'Haar', 
  %   or 'Curvelet'.  Note that curvelet requires the Curvelab library.
  % verbose - if true, prints informative statements
  % waveletType - (deprecated) either 'Daubechies' for Daubechies-4 (default) or 'Haar'
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
  p.addParameter( 'checkAdjoints', false, @islogical );
  p.addParameter( 'debug', false, @(x) isnumeric(x) || islogical(x) );
  p.addParameter( 'nIter', [], @ispositive );
  p.addParameter( 'printEvery', 1, @ispositive );
  p.addParameter( 'transformType', [], @(x) true );
  p.addParameter( 'wavSplit', wavSplit, @isnumeric );
  p.addParameter( 'verbose', false, @(x) isnumeric(x) || islogical(x) );
  p.addParameter( 'waveletType', [], @(x) true );
  p.parse( varargin{:} );
  checkAdjoints = p.Results.checkAdjoints;
  debug = p.Results.debug;
  nIter = p.Results.nIter;
  printEvery = p.Results.printEvery;
  transformType = p.Results.transformType;
  waveletType = p.Results.waveletType;
  verbose = p.Results.verbose;

  if numel( transformType ) == 0
    if numel( waveletType ) > 0
      transformType = waveletType;
    else
      transformType = 'Daubechies';
    end
  end

  if numel( nIter ) == 0
    if debug == true
      nIter = 30;
    else
      nIter = 100;
    end
  end

  M = ( samples ~= 0 );
  b = samples( M == 1 );
  %x0 = zeros( size( samples ) );
  x0 = Fadj( samples );

  % RI = [ Re; Im; ]
  % A = M F RI , A' = (RI)' * F' * M
  % A' * A = (RI)' * F' * M * F * RI
  % gGrad = A'*A*x - A'*b;

  function out = F( x )
    out = fftshift( fftshift( ufft2( ifftshift( ifftshift( x, 1 ), 2 ) ), 1 ), 2 );
  end

  function out = Fadj( y )
    out = fftshift( fftshift( uifft2( ifftshift( ifftshift( y, 1 ), 2 ) ), 1 ), 2 );
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
      cx = fdct_wrapping( x, false );
      out = curvCell2Vec( cx );
    else
      cx = vec2CurvCell( x, curvCellSizes );
      out = ifdct_wrapping( cx, false );
    end
  end

  if strcmp( transformType, 'Curvelet' )
    sparsifier = @(x) curvelet( x );
    sparsifierH = @(x) curvelet( x, 'transp' );

    cx0 = fdct_wrapping( x0, false );
    curvCellSizes = findCellSizes( cx0 );

  elseif strcmp( transformType, 'Haar' )
    sparsifier = @(x) wtHaar2( x, wavSplit );
    sparsifierH = @(y) iwtHaar2( y, wavSplit );

  elseif strcmp( transformType, 'Daubechies' )
    sparsifier = @(x) wtDaubechies2( x, wavSplit );
    sparsifierH = @(y) iwtDaubechies2( y, wavSplit );

  else
    error( 'Unrecognized wavelet type' );
  end
  
  function out = A( x )
    PsiHx = sparsifierH( x );
    FPsiHx = F( PsiHx );
    out = FPsiHx( M == 1 );
  end

  MTy = zeros( size(M) );
  function out = Aadj( y )
    MTy( M == 1 ) = y;
    FadjMTy = Fadj( MTy );
    out = sparsifier( FadjMTy );
  end
  
  if checkAdjoints == true
    x1 = rand( size(samples) ) + 1i * rand( size( samples ) );
    if checkAdjoint( x1, sparsifier, sparsifierH ) ~= true, error( 'sparsifier adjoint is incorrect' ); end
    if checkAdjoint( x1, @F, @Fadj ) ~= true, error( 'Fadj is not the transpose of F' ); end

    xA = sparsifier( x1 );
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

  nPixels = numel( samples );
  proxth = @(x,t) proxL1Complex( x, t * lambda / nPixels );

  function out = h( x )
    out = sum( abs( x(:) ) ) * lambda / nPixels;
  end

  t0 = 1;
  PsiX0 = sparsifier( x0 );
  if debug
    [xStar,oValues] = pogm( PsiX0, @g, @gGrad, proxth, nIter, 'h', @h, 't', t0, ...
      'verbose', verbose, 'printEvery', printEvery );
    %[recon,oValues] = fista( PsiX0, @g, @gGrad, proxth, 'N', nIter, 'h', @h, 'verbose', verbose );   %#ok<ASGLU>
    %[xStar,oValues] = fista_wLS( PsiX0, @g, @gGrad, proxth, 'h', @h, ...
    %  't0', t0, 'N', nIter, 'verbose', true, 'printEvery', printEvery );
  else
    xStar = pogm( PsiX0, @g, @gGrad, proxth, nIter, 't', t0 );
    %recon = fista( PsiX0, @g, @gGrad, proxth, 'N', nIter );   %#ok<UNRCH>
    %xStar = fista_wLS( PsiX0, @g, @gGrad, proxth, 't0', t0, 'N', nIter, ...
    %  'verbose', verbose, 'printEvery', printEvery );
  end

  recon = sparsifierH( xStar );
end

