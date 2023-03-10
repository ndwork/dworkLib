
function [ img, sMaps ] = mri_reconJSENSE( kData, varargin )
  %
  % Inputs:
  % kData - a two dimensional array of complex values; uncollected data have values of 0
  %         Its size is [ nKy nKx nCoils ].
  %
  % Optional Inputs:
  % maxIter - a scalar reprenting the maximum number of iterations
  % polyOrder - either a scalar representing the order in both dimensions or
  %             a two element array representing [ yOrder xOrder ]
  %
  % Outputs:
  % img - a two dimensional complex array that is the reconstructed image
  %
  % Optional Outputs:
  % sMaps - a three dimensional complex array of the sensitivity maps
  %
  % Written by Nicholas Dwork - Copyright 2023
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addParameter( 'maxIter', 10, @(x) ispositive(x) && mod(x,1)==0 );
  p.addParameter( 'polyOrder', [ 17 17 ], @(x) isnonnegative(x) && isinteger(x) );
  p.addParameter( 'relDiffThresh', 0 );
  p.parse( varargin{:} );
  maxIter = p.Results.maxIter;
  polyOrder = p.Results.polyOrder;
  relDiffThresh = p.Results.relDiffThresh;

  if numel( polyOrder ) == 1, polyOrder = [ polyOrder polyOrder ]; end

  img = mri_reconSSQ( kData );

  for iter = 1 : maxIter
    disp([ 'Working on JSENSE iteration ', num2str(iter) ]);

    sMaps = jsense_findSMaps( kData, img, polyOrder );

    img = jsense_findImg( kData, sMaps, img );

    if relDiffThresh > 0
      objDiff = ( objValue - lastObjValue ) / lastObjValue;
      if objDiff < relDiffThresh
        break;
      end
    end
  end
end


function img = jsense_findImg( kData, sMaps, img0 )

  function out = applyS( in, type )
    if nargin < 2 || strcmp( type, 'notransp' )
      out = bsxfun( @times, sMaps, in );
    else
      out = sum( conj(sMaps) .* in, 3 );
    end
  end

  function out = applyF( in, type )
    if nargin < 2 || strcmp( type, 'notransp' )
      out = fftshift2( fft2( ifftshift2( in ) ) );
    else
      out = fftshift2( fft2h( ifftshift2( in ) ) );
    end
  end

  function out = applySF( in, type )
    if nargin < 2 || strcmp( type, 'notransp' )
      out = applyS( applyF( in ) );
    else
      out = applyF( applyS( in, 'transp' ), 'transp' );
    end
  end

  sKData = size( kData );
  dataMask = ( abs(kData) ~= 0 );
  function out = applyE( in, type )
    if nargin < 2 || strcmp( type, 'notransp' )
      in = reshape( in, sKData(1:2) );
      out = applyF( applyS( in ) );
      out = out( dataMask == 1 );
    else
      tmp = zeros( sKData );
      tmp( dataMask == 1 ) = in;
      in = tmp;  clear tmp;
      in = reshape( in, sKData );
      out = applyS( applyF( in .* dataMask, 'transp' ), 'transp' );
    end
    out = out(:);
  end

  doCheckAdjoint = false;
  if doCheckAdjoint == true
    innerProd = @(x,y) real( dotP( x, y ) );
    [checkS,errS] = checkAdjoint( img0, @applyS, 'innerProd', innerProd );   %#ok<ASGLU> 
    [checkF,errF] = checkAdjoint( repmat(img0,[1,1,8]), @applyF, 'innerProd', innerProd );   %#ok<ASGLU> 
    [checkSF,errSF] = checkAdjoint( img0, @applySF, 'innerProd', innerProd );   %#ok<ASGLU> 
    [checkE,errE] = checkAdjoint( img0, @applyE, 'innerProd', innerProd );   %#ok<ASGLU> 
    if checkE ~= 1, error( 'Adjoint check failed' ); end
  end

  img = lsqr( @applyE, kData( dataMask == 1 ), [], 100, [], [], img0(:) );
  img = reshape( img, sKData(1:2) );
end


function sMaps = jsense_findSMaps( kData, img, polyOrder )
  sKData = size( kData );
  nPixels = prod( sKData(1:2) );
  nCoils = sKData( 3 );
  coords = size2imgCoordinates( sKData(1:2) );
  [ ys, xs ] = ndgrid( coords{1} / max(coords{1}), coords{2} / max( coords{2} ) );

  kMask = max( abs( kData ), [], 3 ) ~= 0;

  % Make the 2D Vandermonde matrix
  V = zeros( sKData(1), sKData(2), prod( polyOrder+1 ) );
  orderIndx = 0;
  for xOrder = 0 : polyOrder(2)
    for yOrder = 0 : polyOrder(1)
      orderIndx = orderIndx + 1;
      V(:,:,orderIndx) = ( xs.^xOrder ) .* ( ys.^yOrder );
    end
  end
  Vmat = reshape( V, nPixels, [] );

  function out = applyA( in, type )
    if nargin < 2 || strcmp( type, 'notransp' )
      in = reshape( in, size( V, 3 ), [] );
      Vin = reshape( Vmat * in, sKData );
      imgVin = bsxfun( @times, Vin, img );
      out = fftshift2( fft2( ifftshift2( imgVin ) ) );
      out = reshape( out, nPixels, [] );
      out = out( kMask ~= 0, : );
    else
      kIn = zeros( sKData );
      kIn( abs( kData ) ~= 0 ) = in;
      FhkIn = fftshift2( fft2h( ifftshift2( kIn ) ) );
      imgFhkIn = bsxfun( @times, FhkIn, conj(img) );
      out = Vmat' * reshape( imgFhkIn, nPixels, [] );
    end
    out = out(:);
  end

  doCheckAdjoint = false;
  if doCheckAdjoint == true
    [checkA,errA] = checkAdjoint( rand(324,8), @applyA );
    if checkA ~= true, errro( 'Adjoint check failed' ); end
  end

  polyCoeffs = lsqr( @applyA, kData( abs(kData) ~= 0 ), [], 100, [] );
  polyCoeffs = reshape( polyCoeffs, [], nCoils );

  sMaps = tensorprod( V, polyCoeffs, 3, 1 );
end


