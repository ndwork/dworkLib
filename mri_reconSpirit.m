
function img = mri_reconSpirit( kData, sACR, wSize, varargin )
  % img = mri_reconSpirit( kData, sACR, wSize )
  %
  % Inputs:
  % kData - a complex matrix of size M x N x C, where C is the number of coils
  %   Uncollected data will have a value of 0.
  % sACR - a scalar or a 2 element array that specifies the size of the autocalibration region.
  %  If it's a scalar, then the size is assumed to be [ sACR sACR ].
  % wSize - either a scalar or a two element array of odd vlaues that specifies the size of the interpolation 
  %   kernel.  If wSize is a scalar, then the kernel is assumed to have size [ wSize wSize ].
  %
  % Outputs
  % img - a 2D complex array of size M x N that is the output image.
  %
  % Written by Nicholas Dwork - Copyright 2025
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  if nargin < 3
    disp( 'Usage: img = mri_reconSpirit( kData, sACR, wSize [, ''epsilon'', epsilon ] )' );
    if nargout > 0, img = []; end
    return;
  end

  p = inputParser;
  p.addParameter( 'doChecks', false );
  p.addParameter( 'epsilon', [], @isnonnegative );
  p.addParameter( 'verbose', false );
  p.parse( varargin{:} );
  doChecks = p.Results.doChecks;
  epsilon = p.Results.epsilon;
  verbose = p.Results.verbose;

  if min( mod( wSize, 2 ) ) < 1, error( 'wSize elements must be odd' ); end
  if isscalar( sACR ), sACR = [ sACR sACR ]; end
  if isscalar( wSize ), wSize = [ wSize wSize ]; end
  [ M, N, nCoils ] = size( kData );

  %-- Find the interpolation coefficients w
  acr = cropData( kData, [ sACR(1) sACR(2) nCoils ] );
  w = findW( acr, wSize );

  %-- Use the interpolation coefficients to estimate the missing data
  flipW = padData( flipDims( w, 'dims', [1 2] ), [ M N nCoils nCoils ] );
  function out = applyW( in, op )
    if nargin < 2 || strcmp( op, 'notransp' )
      out = squeeze( sum( circConv2( flipW, in ), 3 ) );
    else
      in = repmat( reshape( in, [ M N 1 nCoils ] ), [ 1 1 nCoils, 1] );
      out = circConv2( flipW, in, 'transp', 'ndimsOut', ndims(kData) );
    end
  end

  if doChecks == true
    [chkW,errChkW] = checkAdjoint( rand( M, N, nCoils ) + 1i * rand( M, N, nCoils ), @applyW );
    if chkW == true
      disp( 'Check of W Adjoint passed' );
    else
      error([ 'Check of W Adjoint failed with error ', num2str(errChkW) ]);
    end
  end

  if numel( epsilon ) == 0 || epsilon == 0
    img = mri_reconSpirit_eps0( @applyW, kData, doChecks );
  else
    img = mri_reconSpirit_epsNonzero( @applyW, kData, epsilon, verbose );
  end
end

function [w, gammas] = findW( acr, wSize )
  %-- Find the interpolation coefficients w
  sACR = size( acr );
  nCoils = sACR( 3 );

  w = cell( 1, 1, 1, nCoils );
  gammas = zeros( nCoils, 1 );
  parfor coilIndx = 1 : nCoils
    A = zeros( (  sACR(2) - wSize(2) + 1 ) * ( sACR(1) - wSize(1) + 1 ), ...
                  wSize(1) * wSize(2) * nCoils - 1 );   %#ok<PFBNS>
    if size( A, 1 ) < size( A, 2 ), error( 'The size of the ACR is too small for this size kernel' ); end
    b = zeros( size(A,1), 1 );
    pt2RemoveIndx = ceil( wSize(1)/2 ) + floor( wSize(2)/2 ) * wSize(1) + ( coilIndx - 1 ) * wSize(2) * wSize(1);
    aIndx = 1;
    for i = ceil( wSize(2)/2 ) : sACR(2) - floor( wSize(2)/2 )
      for j = ceil( wSize(1)/2 ) : sACR(1) - floor( wSize(1)/2 )
        subACR = acr( j - floor( wSize(2)/2 ) : j + floor( wSize(2)/2 ), ...
                               i - floor( wSize(2)/2 ) : i + floor( wSize(2)/2 ), : );   %#ok<PFBNS>
        subACR = subACR(:);
        b( aIndx ) = subACR( pt2RemoveIndx );
        subACR = [ subACR( 1 : pt2RemoveIndx-1 ); subACR( pt2RemoveIndx + 1 : end ); ];
        A( aIndx, : ) = transpose( subACR );
        aIndx = aIndx + 1;
      end
    end
    wCoil = A \ b(:);
    gammas( coilIndx ) = norm( A * wCoil - b )^2 / numel(b);
    wCoil = [ wCoil( 1 : pt2RemoveIndx - 1 ); 0; wCoil( pt2RemoveIndx : end ); ];
    w{1,1,1,coilIndx} = reshape( wCoil, [ wSize nCoils ] );
  end
  w = cell2mat( w );
  %gammas = mean( gammas ) / nCoils;
end

function img = mri_reconSpirit_eps0( applyW, kData, doChecks )
  % minimize || W k - k ||_2 over kEst
  %   where k = toMatrix( D )^T kCollected + ( toMatrix( D^C ) )^T kEst
  %
  %   Here, D is a set of sample indices and k are the sample values that were collected.
  %   D^C is the complement of D
  %   That is, k = kData( D^C ).
  %
  % Equivalently, we want to minimize || ( W - I ) k ||_2 over kEst.
  % Equivalently, minimize || ( W - I ) ( toMatrix( D^C ) )^T kEst + ( W - I ) toMatrix( D )^T kCollected ||_2.
  % Equivalently, minimize || ( W - I ) ( toMatrix( D^C ) )^T kEst + ( W - I ) kData ||_2.

  [ M, N, nCoils ] = size( kData );

  sampleMask = kData ~= 0;
  nSamples = sum( sampleMask(:) );
  nEst = numel( sampleMask ) - nSamples;

  % ( W - I ) ( toMatrix( D^C ) )^T kEst == ( W - I ) kEstMatrix
  kEstMatrix = zeros( M, N, nCoils );
  function out = applyA( in, op )
    if nargin < 2 || strcmp( op, 'notransp' )
      kEstMatrix( not( sampleMask ) ) = in;
      out = applyW( kEstMatrix ) - kEstMatrix;
    else
      in = reshape( in, [ M N nCoils] );
      tmp = applyW( in, 'transp' ) - in;
      out = tmp( not( sampleMask ) );
    end
    out = out(:);
  end

  if doChecks == true
    k0 = rand( nEst, 1 ) + 1i * rand( nEst, 1 );
    [chkA,errChkA] = checkAdjoint( k0, @applyA );
    if chkA == true
      disp( 'Check of A Adjoint passed' );
    else
      error([ 'Check of A Adjoint failed with error ', num2str(errChkA) ]);
    end
  end

  % b = -( W - I ) toMatrix( D )^T kCollected
  b = applyW( kData ) - kData;  b = -b(:);

  k0 = zeros( nEst, 1 );
  tol = 1d-6;
  nMaxIter = 1000;
  [ kStar, lsqrFlag, lsqrRelRes, lsqrIter, lsqrResVec ] = lsqr( @applyA, b, tol, nMaxIter, [], [], k0(:) );   %#ok<ASGLU>

  kOut = kData;
  kOut( not( sampleMask ) ) = kStar;
  img = mri_reconRoemer( mri_reconIFFT( kOut, 'dims', [ 1 2 ] ) );
end


function img = mri_reconSpirit_epsNonzero( applyW, kData, epsilon, verbose )
  % k = kCollected  U  kEst
  %
  % minimize (1/2) || W k - k ||_2^2 over k
  % subject to || D k - kCollected ||_2 <= epsilon
  %   where k = toMatrix( D )^T kCollected + ( toMatrix( D^C ) )^T kEst
  %
  %   Here, D is a set of sample indices and k are the sample values that were collected.
  %   D^C is the complement of D
  %   That is, k = kData( D^C ).
  %
  % Equivalently, we want to minimize || ( W - I ) k ||_2 over kEst.
  % Equivalently, minimize || ( W - I ) ( toMatrix( D^C ) )^T kEst + ( W - I ) toMatrix( D )^T kCollected ||_2.
  % Equivalently, minimize || ( W - I ) ( toMatrix( D^C ) )^T kEst + ( W - I ) kData ||_2.

  function out = g( in )
    out = 0.5 * norm( applyW( in ) - in, 'fro' )^2;
  end

  function out = gGrad( in )
    tmp = applyW( in ) - in;
    out = applyW( tmp, 'transp' ) - tmp;
  end

  sampleMask = kData ~= 0;
  kCollected = kData( sampleMask == 1 );
  function out = h( in )
    diffEst = norm( in( sampleMask == 1 ) - kCollected(:) );
    out = indicatorFunction( diffEst, [ 0, epsilon ] );
  end

  function out = proxth( in, t )   %#ok<INUSD>
    out = in;
    out( sampleMask == 1 ) = projectOntoBall( in( sampleMask == 1 ) - kCollected(:), epsilon ) + kCollected(:);
  end

  [ kStar, objValues ] = fista_wLS( kData, @g, @gGrad, @proxth, 'h', @h, 'verbose', verbose );   %#ok<ASGLU>

  img = mri_reconRoemer( mri_reconIFFT( kStar, 'dims', [ 1 2 ] ) );
end

