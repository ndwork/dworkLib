
function img = mriRecon( kData, varargin )
  % img = mri_Recon( kData [, 'epsData', epsData, 'epsSupport', epsSupport, 'gamma', gamma, 
  %       'lambda', lambda, 'Psi', Psi 'sACR', sACR, 'sMaps', sMaps, 'support', support, 
  %       'verbose', verbose, 'wSize', wSize ] )
  %
  % If nCoils == 1
  %  By default, reconstructs with an inverse Fast Fourier transform
  %
  %  If support is supplied, then solves the following problem
  %    minimize || M F PT x - b ||_2
  %  where P is the matrix that extracts the pixels within the support into a vector
  %
  %  If epsSupport is also supplied, then solves the following problem
  %    minimize || M F x - b ||_2  subject to  Var( Pc x ) <= epsSupport
  %  where Pc is the complement matrix to P.  It extracts the pixels outside the support into a vector.
  %
  %
  % If nCoils > 1
  %   By default, reconstructs with Roemer reconstruction
  %
  % If sMaps is supplied, requires nCoils > 1 and reconstructs by solving the following problem
  %   minimize || M F S x - b ||_2
  %
  % If support is supplied, reconstructs by solving the following problem
  %   minimize || M F PT x - b ||_2
  %
  % If support is supplied, reconstructs by solving the following problem
  %   minimize || M F S PT x - b ||_2
  %
  % Inputs:
  % kData - an array of size M x N x C where C is the number of coils
  %         and uncollected data have value 0.
  %
  % Optional Inputs:
  % sMaps - a complex array of size M x N x C specifying the sensitivity of each coil
  %
  %
  % Outputs
  % img - a 2D complex array of size M x N that is the output image.
  %
  % Written by Nicholas Dwork, Copyright 2025
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addParameter( 'doChecks', true );
  p.addParameter( 'epsData', [], @isnonnegative );
  p.addParameter( 'epsSpiritReg', [], @isnonnegative );
  p.addParameter( 'epsSupport', 0, @isnonnegative );
  p.addParameter( 'lambda', [], @isnonnegative );
  p.addParameter( 'oPsi', true );
  p.addParameter( 'optAlg', [] );
  p.addParameter( 'Psi', [] );
  p.addParameter( 'sACR', [] );
  p.addParameter( 'sMaps', [], @isnumeric );
  p.addParameter( 'support', [] );
  p.addParameter( 'verbose', true );
  p.addParameter( 'wSize', [], @ispositive );
  p.parse( varargin{:} );
  doChecks = p.Results.doChecks;
  epsData = p.Results.epsData;
  epsSpiritReg = p.Results.epsSpiritReg;
  epsSupport = p.Results.epsSupport;
  lambda = p.Results.lambda;
  oPsi = p.Results.oPsi;
  optAlg = p.Results.optAlg;
  Psi = p.Results.Psi;
  sACR = p.Results.sACR;
  sMaps = p.Results.sMaps;
  support = p.Results.support;
  verbose = p.Results.verbose;
  wSize = p.Results.wSize;


  %-- Preliminary variable definitions and input checks
  sampleMask = ( kData ~= 0 );
  sKData = size( kData );
  nKData = prod( sKData );
  sImg = size( kData, [1 2] );
  nImg = prod( sImg );
  nCoils = size( kData, 3 );

  if isscalar( epsData ), epsData = epsData * ones( nCoils, 1 ); end
  if isscalar( sACR ), sACR = [ sACR sACR ]; end
  if isscalar( wSize ), wSize = [ wSize wSize ]; end

  if min( mod( wSize, 2 ) ) < 1, error( 'wSize elements must be odd' ); end
  if numel( wSize ) > 0  &&  numel( sACR ) == 0
    error( 'Must supply sACR if you supply wSize' );
  end
  if numel( wSize ) > 0  &&  nCoils == 1
    error( 'Cannot perform spiritReg with only one coil' );
  end

  if numel( wSize ) > 0
    acr = cropData( kData, [ sACR(1) sACR(2) nCoils ] );
    if numel( epsSpiritReg ) == 0
      [w, epsSpiritReg] = findW( acr, wSize );
    elseif isscalar( epsSpiritReg )
      epsSpiritReg = epsSpiritReg * ones( nCoils, 1 );
    end
    flipW = padData( flipDims( w, 'dims', [1 2] ), [ sImg nCoils nCoils ] );
  end

  if ( numel( lambda ) > 0  &&  lambda > 0 )  ||  numel( epsData ) > 0
    cs = true;
    if ( numel( lambda ) > 0  &&  lambda > 0 )  &&  numel( epsData ) > 0
      error( 'Cannot supply both lambda and epsData')
    end
  else
    cs = false;
  end

  if numel( Psi ) == 0
    wavSplit = makeWavSplit( sImg );
    if oPsi == false
      warning( 'The default value of oPsi being overwritten because default Psi is orthogonal' );
    end
    oPsi = true;
    Psi = @defaultPsi;
  end

  if numel( support ) > 0
    if any( sImg ~= size( support ) ), error( 'Support must be the size of the image' ); end
    nSupport = sum( support(:) );
    nSupportC = numel( support ) - nSupport;
  end

  if numel( sMaps ) > 0  &&  nCoils == 1
    error( 'Cannot supply sMaps with single coil data' );
  end

  coilImgs0 = applyF( kData, 'inv' );
  if nCoils > 1
    img0 = mri_reconRoemer( coilImgs0, 'sMaps', sMaps );
  else
    img0 = coilImgs0;
  end
  sImg0 = size( img0 );

  b = applyM( kData );
  nb = numel( b );
  nbPerCoil = nb / nCoils;

  %-- Function definitions

  function out = applyF( in, op )
    if nargin < 2 || strcmp( op, 'notransp' )
      out = fftshift2( fft2( ifftshift2( in ) ) );
    elseif strcmp( op, 'inv' )
      out = fftshift2( ifft2( ifftshift2( in ) ) );
    elseif strcmp( op, 'transpInv' )
      out = fftshift2( ifft2h( ifftshift2( in ) ) );
    elseif strcmp( op, 'transp' )
      out = fftshift2( fft2h( ifftshift2( in ) ) );
    else
      error( 'unrecognized operation' );
    end
  end

  function out = applyFinv( in, op )
    if nargin < 2  ||  strcmp( op, 'notransp' )
      out = applyF( in, 'inv' );
    else
      out = applyF( in, 'transpInv' );
    end
  end

  function out = applyFS( in, op )
    if nargin < 2 || strcmp( op, 'notransp' )
      out = applyF( applyS( in ) );
    else
      out = applyS( applyF( in, 'transp' ), 'transp' );
    end
  end

  function out = applyM( in, op )
    if nargin < 2  ||  strcmp( op, 'notransp' )
      out = in( sampleMask == 1 );
    else
      out = zeros( sKData );
      out( sampleMask == 1 ) = in;
    end
  end

  function out = applyMF( in, op )
    if nargin < 2  ||  strcmp( op, 'notransp' )
      out = applyM( applyF( in ) );
    else
      out = applyF( applyM( in, op ), op );
    end
  end

  function out = applyMFPT( in, op )
    if nargin < 2  ||  strcmp( op, 'notransp' )
      out = applyMF( applyP( in, 'transp' ) );
    else
      out = applyP( applyMF( in, op ), 'notransp' );
    end
  end

  function out = applyMFPcoilsT( in, op )
    if nargin < 2  ||  strcmp( op, 'notransp' )
      out = applyMF( applyPcoils( in, 'transp' ) );
    else
      out = applyPcoils( applyMF( in, op ), 'notransp' );
    end
  end

  function out = applyMFS( in, op )
    if nargin < 2  ||  strcmp( op, 'notransp' )
      out = applyMF( applyS( in ) );
    else
      out = applyS( applyMF( in, op ), op );
    end
  end

  function out = applyMFSPT( in, op )
    if nargin < 2  ||  strcmp( op, 'notransp' )
      out = applyMFS( applyP( in, 'transp' ) );
    else
      out = applyP( applyMFS( in, op ), 'notransp' );
    end
  end

  function out = applyP( in, op )
    if nargin < 2  ||  strcmp( op, 'notransp' )
      out = in( support == 1 );
    else
      out = zeros( sImg );
      out( support == 1 ) = in;
    end
  end

  function out = applyPcoils( in, op )
    if nargin < 2  ||  strcmp( op, 'notransp' )
      out = zeros( nSupport, nCoils );
      for c = 1 : nCoils
        out(:,c) = applyP( in(:,:,c) );
      end
    else
      out = zeros( sKData );
      for c = 1 : nCoils
        out(:,:,c) = applyP( in(:,c), op );
      end
    end
  end

  function out = applyPC( in, op )
    if nargin < 2  ||  strcmp( op, 'notransp' )
      out = in( support == 0 );
    else
      out = zeros( sImg );
      out( support == 0 ) = in;
    end
  end

  function out = applyPCcoils( in, op )
    if nargin < 2  ||  strcmp( op, 'notransp' )
      out = zeros( nSupportC, nCoils );
      for c = 1 : nCoils
        out(:,c) = applyPC( in(:,:,c) );
      end
    else
      out = zeros( sKData );
      for c = 1 : nCoils
        out(:,:,c) = applyPC( in(:,c), op );
      end
    end
  end

  function out = applyS( in, op )
    if nargin < 2 || strcmp( op, 'notransp' )
      out = bsxfun( @times, sMaps, in );
    else
      out = sum( conj(sMaps) .* in, 3 );
    end
  end

  function out = applyW( in, op )
    if nargin < 2 || strcmp( op, 'notransp' )
      out = squeeze( sum( circConv2( flipW, in ), 3 ) );
    else
      in = repmat( reshape( in, [ sImg 1 nCoils ] ), [ 1 1 nCoils, 1] );
      out = circConv2( flipW, in, 'transp', 'ndimsOut', ndims(kData) );
    end
  end

  function out = applyWmI( in, op )
    % apply W minus I
    if nargin < 2 || strcmp( op, 'notransp' )
      out = applyW( in ) - in;
    else
      out = applyW( in, 'transp' ) - in;
    end
  end

  function out = applyWmIFS( in, op )
    if nargin < 2 || strcmp( op, 'notransp' )
      FSin = applyFS( in );
      out = applyW( FSin ) - FSin;
    else
      out = applyFS( applyW( in, 'transp' ) - in, 'transp' );
    end
  end

  function out = defaultPsi( in, op )
    out = zeros( size( in ) );
    if nargin < 2 || strcmp( op, 'notransp' )
      for sliceIndx = 1 : size( in, 3 )
        out(:,:,sliceIndx) = wtDaubechies2( in(:,:,sliceIndx), wavSplit );
      end
    else
      for sliceIndx = 1 : size( in, 3 )
        out(:,:,sliceIndx) = iwtDaubechies2( in(:,:,sliceIndx), wavSplit );
      end
    end
  end

  function out = indicatorEpsData( in )
    out = 0;
    bTmp = reshape( b, [], nCoils );
    in = reshape( in, [], nCoils );
    for c = 1 : nCoils
      out = indicatorFunction( norm( in - bTmp(:,c) )^2 / ( nbPerCoil-1 ), [ 0 epsData(c) ] );
      if out ~= 0, break; end
    end
  end

  epsBallRadiuses = sqrt( epsData * (nbPerCoil-1) );
  function out = projectOntoEpsBalls( in )
    % proximal operator of indicator( || x ||_2 /  <= ballRadius(c) )
    in = reshape( in, nbPerCoil, nCoils );
    out = zeros( nbPerCoil, nCoils );
    for c = 1 : nCoils
      out(:,c) = projectOntoBall( in(:,c), epsBallRadiuses(c) );
    end
    out = out(:);
  end

  if numel( support ) > 0  &&  numel( epsSupport ) > 0
    supportBallRadius = sqrt( epsSupport * (nSupportC - 1) );
  end
  function out = projOutsideSupportOntoBall( in )
    out = in;
    out( support == 0 ) = projectOntoBall( applyPC( in ), supportBallRadius );
  end

  if numel( wSize ) > 0
    ballRadiusesSpiritReg = sqrt( epsSpiritReg * ( nImg - 1 ) );
  end
  function out = proxSpiritReg( in, t )   %#ok<INUSD>
    % indicator( || in ||_2^2 / (nImg-1) <= epsSpiritReg(coilIndx) ) is equivalent to
    % indicator( || in(:,:,coilIndx) ||_Fro <= sqrt( epsSpiritReg( coilIndx ) * (nImg-1) ) )
    in = reshape( in, sKData );
    out = zeros( sKData );
    for c = 1 : nCoils
      out(:,:,c) = projectOntoBall( in(:,:,c), ballRadiusesSpiritReg(c), 'fro' );
    end
    out = out(:);
  end

  function out = spiritReg( x )
    % indicator( || x ||_2^2 / (nImg-1) <= epsSpiritReg(coilIndx) ) is equivalent to
    % indicator( || x(:,:,coilIndx) ||_Fro <= sqrt( epsSpiritReg( coilIndx ) * (nImg-1) ) )
    x = reshape( x, sKData );
    xNorms = LpNorms( reshape( x, [], nCoils ), 2, 1 );
    out = 0;
    if any( ( xNorms.^2 / ( nImg - 1 ) ) > ballRadiusesSpiritReg )
      out = Inf;
    end
  end

  function out = metricSpiritRegViolation( in )
    WmIFSin = applyWmIFS(in);
    violations = zeros(nCoils,1);
    for c = 1 : nCoils
      violations(c) = abs( norm( WmIFSin(:,:,c), 'fro' ) - ballRadiusesSpiritReg(c) );
    end
    out = max( violations );
  end


  %-- Concatenated functions

  function out = concat_MFS_Pc( in, op )
    if nargin < 2  ||  strcmp( op, 'notransp' )
      out = [ reshape( applyMFS(in), [], 1 );  applyPC(in) ];
    else
      out = applyMFS( in(1:nb), op ) + applyPC( in(nb+1:end), op );
    end
  end

  function out = concat_MFS_WmIFS( in, op )
    if nargin < 2 || strcmp( op, 'notransp' )
      FSin = applyFS( in );
      MFSin = applyM( FSin );
      WmIFSin = applyWmI( FSin );
      out = [ MFSin(:); WmIFSin(:); ];
    else
      n1 = nb;
      n2 = n1 + nKData;
      MTin1 = reshape( applyM( in( 1 : n1 ), op ), sKData );
      WHmIHin3 = reshape( applyWmI( reshape( in( n1+1 : n2 ), sKData ), op ), sKData );
      out = applyFS( MTin1 + WHmIHin3, op );
    end
  end

  function out = concat_MFS_WmIFS_I( in, op )
    if nargin < 2 || strcmp( op, 'notransp' )
      FSin = applyFS( in );
      MFSin = applyM( FSin );
      WmIFSin = applyWmI( FSin );
      out = [ MFSin(:); WmIFSin(:); in(:); ];
    else
      n1 = nb;
      n2 = n1 + nKData;
      n3 = n2 + nImg;
      MTin1 = reshape( applyM( in( 1 : n1 ), op ), sKData );
      WHmIHin2 = reshape( applyWmI( reshape( in( n1+1 : n2 ), sKData ), op ), sKData );
      in3 = reshape( in( n2+1 : n3 ), sImg );
      out = applyFS( MTin1 + WHmIHin2, op ) + in3;
    end
  end

  function out = concat_Psi_MFS_Pc( in, op )
    if nargin < 2  ||  strcmp( op, 'notransp' )
      out = [ reshape( Psi(in), [], 1 ); ...
              reshape( applyMFS(in), [], 1 );
              applyPC( in ) ];
    else
      n1 = nPsi;
      n2 = nPsi + nb;
      n3 = nPsi + nb + ( nImg - nSupport );
      out = Psi( in(1:nPsi), op ) + applyMFS( in(n1:n2), op ) + applyPC( in(n2:n3), op );
    end
  end

  function out = concat_Psi_MFS_WmIFS( in, op )
    if nargin < 2 || strcmp( op, 'notransp' )
      PsiIn = Psi( in );
      FSin = applyFS( in );
      MFSin = applyM( FSin );
      WmIFSin = applyWmI( FSin );
      out = [ PsiIn(:); MFSin(:); WmIFSin(:); ];
    else
      n0 = nPsi;
      n1 = n0 + nb;
      n2 = n1 + nKData;
      PsiIn0 = Psi( reshape( in( 1 : n0 ), sKData), op );
      MTin1 = reshape( applyM( in( n0 : n1 ), op ), sKData );
      WHmIHin3 = reshape( applyWmI( reshape( in( n1+1 : n2 ), sKData ), op ), sKData );
      out = applyFS( PsiIn0 + MTin1 + WHmIHin3, op );
    end
  end


  %-- Vectorized functions

  function out = apply_vMF( in, op )
    if nargin < 2  ||  strcmp( op, 'notransp' )
      in = reshape( in, sKData );
      out = applyMF( in );
    else
      out = applyMF( in, op );
      out = out(:);
    end
  end

  function out = apply_vMFPcoilsT( in, op )
    if nargin < 2  ||  strcmp( op, 'notransp' )
      in = reshape( in, [ nSupport, nCoils ] );
      out = applyMFPcoilsT( in );
    else
      out = applyMFPcoilsT( in, op );
      out = out(:);
    end
  end

  function out = apply_vMFPT( in, op )
    if nargin < 2  ||  strcmp( op, 'notransp' )
      in = reshape( in, [ nSupport, nCoils ] );
      out = applyMFPT( in );
    else
      out = applyMFPT( in, op );
    end
  end

  function out = apply_vMFS( in, op )
    if nargin < 2  ||  strcmp( op, 'notransp' )
      in = reshape( in, sImg );
      out = applyMFS( in );
    else
      out = applyMFS( in, op );
      out = out(:);
    end
  end

  function out = apply_vMFSPT( in, op )
    if nargin < 2  ||  strcmp( op, 'notransp' )
      in = reshape( in, [ nSupport, 1 ] );
      out = applyMFSPT( in );
    else
      out = applyMFSPT( in, op );
      out = out(:);
    end
  end

  function out = vPsi( in, op )
    if numel( in ) == nImg
      in = reshape( in, sImg );
    else
      in = reshape( in, sKData );
    end
    if nargin < 2 || strcmp( op, 'notransp' )
      out = Psi( in );
    else
      out = Psi( in, op );
    end
    out = out(:);
  end


  %-- Checks
  if doChecks
    [chk_M,errChk_M] = checkAdjoint( rand( sKData ) + 1i * rand( sKData ), @applyM );
    if chk_M == true
      disp( 'Check of applyM passed' );
    else
      error([ 'Check of applyM Adjoint failed with error ', num2str(errChk_M) ]);
    end

    [chk_F,errChk_F] = checkAdjoint( rand( sKData ) + 1i * rand( sKData ), @applyF );
    if chk_F == true
      disp( 'Check of applyF passed' );
    else
      error([ 'Check of applyF Adjoint failed with error ', num2str(errChk_F) ]);
    end

    [chk_Finv,errChk_Finv] = checkAdjoint( rand( sKData ) + 1i * rand( sKData ), @applyFinv );
    if chk_Finv == true
      disp( 'Check of applyFinv passed' );
    else
      error([ 'Check of applyFinv Adjoint failed with error ', num2str(errChk_Finv) ]);
    end

    [chk_MF,errChk_MF] = checkAdjoint( rand( sKData ) + 1i * rand( sKData ), @applyMF );
    if chk_MF == true
      disp( 'Check of applyMF passed' );
    else
      error([ 'Check of applyMF Adjoint failed with error ', num2str(errChk_MF) ]);
    end

    [chk_vPsi,errChk_vPsi] = checkAdjoint( rand( sKData ) + 1i * rand( sKData ), @vPsi );
    if chk_vPsi == true
      disp( 'Check of vPsi Adjoint passed on 3D data' );
    else
      error([ 'Check of vPsi Adjoint failed on 3D data with error ', num2str(errChk_vPsi) ]);
    end

    [chk_vPsi,errChk_vPsi] = checkAdjoint( rand( sImg ) + 1i * rand( sImg ), @vPsi );
    if chk_vPsi == true
      disp( 'Check (2nd) of vPsi Adjoint passed on 2D data' );
    else
      error([ 'Check (2nd) of vPsi Adjoint failed on 2D data with error ', num2str(errChk_vPsi) ]);
    end

    if numel( support ) > 0
      [chkP,errChkP] = checkAdjoint( rand( sImg ) + 1i * rand( sImg ), @applyP );
      if chkP == true
        disp( 'Check of P adjoint passed' );
      else
        error([ 'Check of P adjoint failed with error ', num2str(errChkP) ]);
      end
  
      [chkPc,errChkPc] = checkAdjoint( rand( sImg ) + 1i * rand( sImg ), @applyPC );
      if chkPc == true
        disp( 'Check of Pc adjoint passed' );
      else
        error([ 'Check of Pc adjoint failed with error ', num2str(errChkPc) ]);
      end

      if nCoils == 1

        [chk_vMFPT,errChk_vMFPT] = checkAdjoint( rand(nSupport,nCoils) + 1i * rand(nSupport,nCoils), ...
          @apply_vMFPT );
        if chk_vMFPT == true
          disp( 'Check of vMFPT adjoint passed' );
        else
          error([ 'Check of vMFPT adjoint failed with error ', num2str(errChk_vMFPT) ]);
        end

      else

        [chkPcoils,errChkPcoils] = checkAdjoint( rand( sKData ) + 1i * rand( sKData ), @applyPcoils );
        if chkPcoils == true
          disp( 'Check of Pcoils adjoint passed' );
        else
          error([ 'Check of Pcoils adjoint failed with error ', num2str(errChkPcoils) ]);
        end

        [chkPCcoils,errChkPCcoils] = checkAdjoint( rand( sKData ) + 1i * rand( sKData ), @applyPCcoils );
        if chkPCcoils == true
          disp( 'Check of PCcoils adjoint passed' );
        else
          error([ 'Check of PCcoils adjoint failed with error ', num2str(errChkPCcoils) ]);
        end
      end
    end

    if numel( wSize ) > 0
      [chkWmI,errChkWmI] = checkAdjoint( rand( sKData ) + 1i * rand( sKData ), @applyWmI );
      if chkWmI == true
        disp( 'Check of W minus I Adjoint passed' );
      else
        error([ 'Check of W minus I Adjoint failed with error ', num2str(errChkWmI) ]);
      end
    end

    if numel( sMaps ) > 0
      [chk_vMFS,errChk_vMFS] = checkAdjoint( rand( sImg ) + 1i * rand( sImg ), @apply_vMFS );
      if chk_vMFS == true
        disp( 'Check of vMFS adjoint passed' );
      else
        error([ 'Check of vMFS adjoint failed with error ', num2str(errChk_vMFS) ]);
      end

      if numel( wSize ) > 0
        [chkWmIFS,errChkWmIFS] = checkAdjoint( rand( sImg0 ) + 1i * rand( sImg0 ), @applyWmIFS );
        if chkWmIFS == true
          disp( 'Check of WmIFS Adjoint passed' );
        else
          error([ 'Check of WmIFS Adjoint failed with error ', num2str(errChkWmIFS) ]);
        end
      end

      if numel( support ) > 0
        tmp = rand( nSupport, 1 ) + 1i * rand( nSupport, 1 );
        [chk_vMFSPT,errChk_vMFSPT] = checkAdjoint( tmp(:), @apply_vMFSPT );
        if chk_vMFSPT == true
          disp( 'Check of vMFSPT adjoint passed' );
        else
          error([ 'Check of vMFSPT adjoint failed with error ', num2str(errChk_vMFSPT) ]);
        end

        [ chk_MFS_Pc, errChk_MFS_Pc ] = checkAdjoint( rand( sImg ) + 1i * rand( sImg ), @concat_MFS_Pc );
        if chk_MFS_Pc == true
          disp( 'Check of concat_MFS_Pc adjoint passed' );
        else
          error([ 'Check of concat_MFS_Pc adjoint failed with error ', num2str(errChk_MFS_Pc) ]);
        end
      end
    end
  end  % if doChecks

  %-- Reconstruct the image

  tol = 1d-6;
  nMaxIter = 1000;

  if cs == true

    nPsiCoils = numel( Psi( coilImgs0 ) );
    nPsi = numel( Psi( img0 ) );

    if nCoils > 1
  
      if numel( sMaps ) > 0

        if numel( support ) > 0

          if epsSupport > 0

            if numel( wSize ) > 0
              % minimize || Psi x ||_1
              % subject to  || M F S(c) x - b ||_2^2 / (nbPerCoil-1) <= epsData(c)
              %        and  || Pc x ||_2^2 / ( nImg - Ns - 1 ) <= epsSupport
              %        and  || ((W-I) F S x)(:,:,c) ||_Fro <= sqrt( epsSpiritReg(c) * nImg ) for all c

              g1 = @( in ) indicatorEpsData( in );
              proxg1 = @(in,t) projectOntoEpsBalls( in - b ) + b;

              if oPsi == true
                applyA = @(in,op) concat_MFS_WmIFS_I( in, op );
                f = @(in) norm( reshape( Psi( in ), [], 1 ), 1 );
                proxf = @(in,t) proxCompositionAffine( @proxL1Complex, in, Psi, 0, 1, t );
                g2 = @spiritReg;
                proxg2 = @proxSpiritReg;
                g3 = @(in) indicatorFunction( ...
                  norm( in, 2 )^2 / ( nImg - nSupport - 1 ), [ 0 epsSupport ] );
                proxg3 = @(in,t) projOutsideSupportOntoBall( in );
                n1 = nb;
                n2 = n1 + nKData;
                n3 = n2 + nImg;
                g = @(in) g1( in(1:n1) ) + g2( in(n1+1:n2) ) + g3( in(n2+1:n3) );
                proxg = @(in,t) [ proxg1( in(1:n1), t );  ...
                                  proxg2( in(n1+1:n2) ); ...
                                  proxg3( in(n2+1:n3) ); ];
                proxgConj = @(in,s) proxConj( proxg, in, s );
    
                metric1 = @(in) 0.5 * norm( applyMFS( in ) - b )^2;
                metric2 = @metricSpiritRegViolation;
                metrics = { metric1, metric2 };
                metricNames = { 'sparsity', 'spirit violation' };
                [img, objValues, mValues] = pdhgWLS( img0, proxf, proxgConj, 'A', applyA, 'f', f, 'g', g, ...
                  'N', nMaxIter, 'metrics', metrics, 'metricNames', metricNames, 'printEvery', 10, 'verbose', verbose );   %#ok<ASGLU>

              else
                applyA = @(in,op) concat_Psi_MFS_WmIFS( in, op );
                f = @(in) indicatorFunction( ...
                  norm( in, 2 )^2 / ( nImg - nSupport - 1 ), [ 0 epsSupport ] );
                proxf = @(in,t) projOutsideSupportOntoBall( in );
                g0 = @(in) norm( in(:), 1 );
                proxg0 = @proxL1Complex;
                g2 = @spiritReg;
                proxg2 = @proxSpiritReg;
                n0 = nKData;
                n1 = n0 + nb;
                n2 = n1 + nKData;
                g = @(in) g0( in(1:n0) ) + g1( in(n0:n1) ) + g2( in(n1+1:n2) );
                proxg = @(in,t) [ proxg0( in(1:n0), t ); ...
                                  proxg1( in(n0+1:n1), t ); ...
                                  proxg2( in(n1+1:n2) ) ];
                proxgConj = @(in,s) proxConj( proxg, in, s );

                metric1 = @(in) 0.5 * norm( applyMFS( in ) - b )^2;
                metric2 = @metricSpiritRegViolation;
                metrics = { metric1, metric2 };
                metricNames = { 'sparsity', 'spirit violation' };
                [img, objValues, mValues] = pdhgWLS( img0, proxf, proxgConj, 'A', applyA, 'f', f, 'g', g, ...
                  'N', nMaxIter, 'metrics', metrics, 'metricNames', metricNames, 'printEvery', 10, 'verbose', verbose );   %#ok<ASGLU>
              end
            else

              if numel( epsData ) > 0
                % minimize || Psi x ||_1  subject to  || M F S(c) x - b(c) ||_2^2 / ( nbPerCoil - 1 ) <= epsData(c)
                %                                and  || Pc x ||_2^2 / ( nImg - Ns - 1 ) <= epsSupport
  
                g1 = @( in ) indicatorEpsData( in );
                proxg1 = @(in,t) projectOntoEpsBalls( in - b ) + b;

                if oPsi == true
                  f = @(in) normL2L1( Psi( in ) );
                  proxf = @(in,t) proxCompositionAffine( @proxL2L1, in, Psi, 0, 1, t );
                  applyA = @concat_MFS_Pc;
                  g2 = @(in) indicatorFunction( ...
                    norm( in, 2 )^2 / ( nImg - nSupport - 1 ), [ 0 epsSupport ] );
                  g = @(in) g1( in(1:nb) ) + g2( in(nb+1:end) );
                  proxg2 = @(in,t) projOutsideSupportOntoBall( in, sqrt( epsSupport * (nImg-nSupport-1) ) );
                  proxg = @(in,t) [ proxg1( in(1:nb), t );  proxg2( in(nb+1:end), t ); ];
                  proxgConj = @(in,s) proxConj( proxg, in, s );
                  [img, objValues, mValues] = pdhgWLS( img0, proxf, proxgConj, 'A', applyA, ...
                    'f', f, 'g', g, 'verbose', verbose );   %#ok<ASGLU>
  
                else  % if oPsi == true
                  applyA = @concat_Psi_MFS_Pc;
                  f = @(in) indicatorFunction( ...
                    norm( in(nb+1:end), 2 )^2 / ( nImg - nSupport - 1 ), [ 0 epsSupport ] );
                  proxf = @(in,t) projOutsideSupportOntoBall( in );
                  g0 = @(in) normL2L1( in );
                  n1 = nImg;
                  n2 = n1 + nb;
                  g = @(in) g0(1:n1) + g1( in(n1+1:n2) );
                  proxg0 = @proxL2L1;
                  proxg = @(in,t) [ proxg0( in(1:n1) ); ...
                                    proxg1( in(n1+1:n2), t ); ];
                  proxgConj = @(in,s) proxConj( proxg, in, s );
                  [img, objValues, mValues] = pdhgWLS( img0, proxf, proxgConj, 'A', applyA, ...
                    'f', f, 'g', g, 'verbose', verbose );   %#ok<ASGLU>
                end
  
              else  % if numel( epsData ) > 0
                % lambda provided
                % minimize (1/2) || M F S x - b ||_2^2 + lambda || Psi x ||_1
                % subject to || Pc x ||_2^2 / ( Nc - 1 ) <= epsSupport
                % Could use PD3O based on whether or not oPsi is true
  
                error( 'Not yet implemented' );
  
              end
            end

          else  % if epsSupport > 0

            error( 'Not yet implemented' );
          end

        else  % if numel( support ) > 0

          if numel( wSize ) > 0

            if numel( lambda ) > 0  &&  lambda > 0 
              % minimize (1/2) || M F S x - b ||_2^2  + lambda || Psi ||_1
              % subject to  || ((W-I) F S x)(:,:,c) ||_Fro <= sqrt( epsSpiritReg(c) * nImg ) for all c
              error( 'Not yet implemented' );

            else
              % minimize || Psi x ||_1
              % subject to  || M F S(c) x - b ||_2^2 / (nbPerCoil-1) <= epsData(c)
              %        and  || ((W-I) F S x)(:,:,c) ||_Fro <= sqrt( epsSpiritReg(c) * nImg ) for all c

              g1 = @( in ) indicatorEpsData( in );
              proxg1 = @(in,t) projectOntoEpsBalls( in - b ) + b;
              g2 = @spiritReg;
              proxg2 = @proxSpiritReg;

              metric1 = @(in) 0.5 * norm( applyMFS( in ) - b )^2;
              metric2 = @metricSpiritRegViolation;
              metrics = { metric1, metric2 };
              metricNames = { 'sparsity', 'violation' };

              if oPsi == true
                applyA = @(in,op) concat_MFS_WmIFS( in, op );
                f = @(in) norm( reshape( Psi( in ), [], 1 ), 1 );
                proxf = @(in,t) proxCompositionAffine( @proxL1Complex, in, Psi, 0, 1, t );
                % || M F S x - b ||_2^2 / (nb-1) <= epsData  <==>  || M F S x - b ||_2 <= sqrt( epsData * (nb-1) )
                n1 = nb;
                n2 = n1 + nKData;
                g = @(in) g1( in(1:n1) ) + g2( in(n1+1:n2) );
                proxg = @(in,t) [ proxg1( in(1:n1), t );  proxg2(in(n1+1:n2)) ];
                proxgConj = @(in,s) proxConj( proxg, in, s );
    
                [img, objValues, mValues] = pdhgWLS( img0, proxf, proxgConj, 'A', applyA, 'f', f, 'g', g, ...
                  'N', 1000, 'metrics', metrics, 'metricNames', metricNames, 'printEvery', 10, 'verbose', verbose );   %#ok<ASGLU>

              else
                applyA = @(in,op) concat_Psi_MFS_WmIFS( in, op );
                f = [];
                proxf = [];
                g0 = @(in) norm( in, 1 );
                proxg0 = @proxL1Complex;
                n0 = nKData;
                n1 = n0 + nb;
                n2 = n1 + nKData;
                g = @(in) g0( in(1:n0) ) + g1( in(n0:n1) ) + g2( in(n1+1:n2) );
                proxg = @(in,t) [ proxg0( in(n0:n1), t ); ...
                                  proxg1( in(1:n1), t ); ...
                                  proxg2(in(n1+1:n2)) ];
                proxgConj = @(in,s) proxConj( proxg, in, s );

                [img, objValues, mValues] = pdhgWLS( img0, proxf, proxgConj, 'A', applyA, 'f', f, 'g', g, ...
                  'N', 1000, 'metrics', metrics, 'metricNames', metricNames, 'printEvery', 10, 'verbose', verbose );   %#ok<ASGLU>
              end
            end

          else % if numel( wSize ) > 0

            if oPsi == true
  
              if numel( lambda ) > 0  &&  lambda > 0
                % minimize (1/2) || M F S x - b ||_2^2 + lambda || Psi x ||_1
    
                applyA = @applyMFS;
                f = @(in) 0.5 * norm( applyA(in) - b, 2 )^2;
                ATb = applyA( b, 'transp' );
                fGrad = @(in) applyA( applyA( in ), 'transp' ) - ATb;
                g = @(in) lambda * norm( vPsi( in ), 1 );
                proxg = @(in,t) proxCompositionAffine( @proxL1Complex, in, Psi, 0, 1, lambda * t );
    
                ATA = @(in) applyA( applyA( in ), 'transp' );
                L = powerIteration( ATA, rand( size( img0 ) ), 'symmetric', true );
                t0 = 1 / L;
                if verbose == true
                  [ img, objValues, relDiffs] = fista_wLS( img0, f, fGrad, proxg, 'h', g, 't0', t0, ...
                    'verbose', verbose );   %#ok<ASGLU>
                else
                  img = fista_wLS( img0, f, fGrad, proxg, 'h', g, 'verbose', verbose );
                end
  
              else  % if numel( lambda ) > 0  &&  lambda > 0
                % minimize || Psi x ||_1  subject to  || M F S(c) x - b ||_2^2 / ( nbPerCoil - 1 ) <= epsData(c)

                applyA = @applyMFS;
                f = @(in) norm( vPsi( in ), 1 );
                proxf = @(in,t) proxCompositionAffine( @proxL1Complex, in, Psi, [], 1, t );
                g = @( in ) indicatorEpsData( in );
                proxg = @(in,t) projectOntoEpsBalls( in - b ) + b;
                proxgConj = @(in,s) proxConj( proxg, in, s );
  
                metric1 = f;
                metric2 = @(in) abs( norm( applyA( in ) - b, 2 ) - ballRadius );
                metrics = { metric1, metric2 };
                metricNames = { 'sparsity', 'violation' };
                [img, objValues, mValues] = pdhgWLS( img0, proxf, proxgConj, 'A', applyA, 'f', f, 'g', g, ...
                  'metrics', metrics, 'metricNames', metricNames, 'printEvery', 10, 'verbose', verbose );   %#ok<ASGLU>
  
              end

            else  % if oPsi == true
              if numel( lambda ) > 0  &&  lambda > 0
              else  % if numel( lambda ) > 0  &&  lambda > 0
              end
              error( 'Not yet implemented' );
            end
          end
        end

      else  % if numel( sMaps ) > 0
  
        if numel( support ) > 0

          if epsSupport > 0

          else

          end

        else % if numel( support ) > 0

          if numel( epsData ) > 0
            % minimize || Psi x ||_{2,1} subject to || M F x - b ||_2^2 / ( npPerCoil - 1 ) < epsData(c)

            if oPsi == true
              applyA = @applyMF;
              f = @(in) normL2L1( Psi( in ) );
              proxf = @(in,t) proxCompositionAffine( @proxL2L1, in, Psi, 0, 1, t );
              g = @( in ) indicatorEpsData( in );
              proxg = @(in,t) projectOntoEpsBalls( in - b ) + b;
              proxgConj = @(in,s) proxConj( proxg, in, s );
              [y, objValues, mValues] = pdhgWLS( img0, proxf, proxgConj, 'A', applyA, ...
                'beta', 1, 'tau', tau, 'f', f, 'g', g, 'N', N, 'dsc', true, ...
                'printEvery', 10, 'metrics', metrics, 'metricNames', metricNames, 'verbose', verbose );   %#ok<ASGLU>

            else  % if oPsi == true
              applyA = @concat_Psi_MF;
              f = [];
              proxf = [];
              g1 = @(in) normL2L1( in );
              g2 = @( in ) indicatorEpsData( in );
              g = @(in) g1( reshape(in(1:nKData), sKData) ) + g2( in(nKData+1:end) );
              proxg1 = @(in,t) proxL2L1( reshape( in, sKData ), t );
              proxg2 = @(in,t) projectOntoEpsBalls( in - b ) + b;
              proxg = @(in,t) [ proxg1( in(1:nKData), t ); ...
                                proxg2( in(nKData+1:end), t ); ];
              proxgConj = @(in,s) proxConj( proxg, in, s );
              [y, objValues, mValues] = pdhgWLS( img0, proxf, proxgConj, 'A', applyA, ...
                'beta', 1, 'tau', tau, 'f', f, 'g', g, 'N', N, 'dsc', true, ...
                'printEvery', 10, 'metrics', metrics, 'metricNames', metricNames, 'verbose', verbose );   %#ok<ASGLU>

            end

          else

          end

        end

      end

    else  % if nCoils > 1

      if numel( support ) > 0

        if epsSupport > 0

        else  % if epsSupport > 0

        end

      else  % if numel( support ) > 0

        if numel( lambda ) > 0  &&  lambda > 0
          % minimize (1/2) || M F x - b ||_2^2 + lambda || Psi x ||_1

          if oPsi == true
            g = @(in) 0.5 * norm( applyMF( in ) - b, 2 )^2;
            FTMTb = applyF( kData, 'transp' );
            gGrad = @(in) applyF( applyF( in ) .* sampleMask, 'transp' ) - FTMTb;
            proxth = @(in,t) proxCompositionAffine( in, Psi, [], 1, t * lambda );
            [img,objectiveValues,relDiffs] = fista_wLS( img0, g, gGrad, proxth );   %#ok<ASGLU>

          else
          end

        else
          % minimize || Psi x ||_1  subject to  || M F x - b ||_2^2 / (nbPerCoil - 1) <= epsData(c)
          % if A = M F the A AT = M F FT MT = a scaled identity matrix
          if oPsi == true
            % Solve this problem with ADMM

          else
          end
        end

      end

    end

  else  % if cs == true

    if nCoils > 1
  
      if numel( sMaps ) > 0
  
        if numel( support ) > 0
  
          if epsSupport > 0
            % minimize || M F S x - b ||_2 subject to || Pc x ||_2^2 / ( Nc - 1 ) <= epsSupport
            g = @(in) 0.5 * norm( applyMFS( in ) - b, 2 )^2;
            STFTMTb = applyFS( kData, 'transp' );
            %gGrad = @(in) applyMFS( applyMFS( in ), 'transp' ) - STFTMTb;
            gGrad = @(in) applyFS( applyFS( in ) .* sampleMask, 'transp' ) - STFTMTb;
            proxth = @(in) projOutsideSupportOntoBall( in, sqrt( epsSupport * (nSupportC - 1) ) );
            [img,objectiveValues,relDiffs] = fista_wLS( img0, g, gGrad, proxth );   %#ok<ASGLU>
  
          else  % if epsSupport > 0

            % minimize || M F S PT x - b ||_2
            img0 = applyP( img0 );
            b = applyM( kData );
            [ img, lsqrFlag, lsqrRelRes, lsqrIter, lsqrResVec ] = lsqr( ...
              @apply_vMFSPT, b, tol, nMaxIter, [], [], img0(:) );   %#ok<ASGLU>
            img = applyP( img, 'transp' );
  
          end
  
        else  % if numel( support ) > 0
  
          if numel( wSize ) > 0  % do SPIRiT Regularization
            % minimize (1/2) || M F S x - b ||_2^2  
            % subject to  || ((W-I) F S x)(:,:,c) ||_Fro <= sqrt( epsSpiritReg(c) * nImg ) for all c
            applyA = @(in,op) concat_MFS_WmIFS( in, op );
            n1 = nb;
            n2 = n1 + nKData;
            f = [];
            proxf = [];
            g1 = @( x ) 0.5 * norm( x - b )^2;
            proxg1 = @(in,t) proxL2Sq( in, t, b );
            g2 = @spiritReg;
            proxg2 = @proxSpiritReg;
            g = @(in) g1( in(1:n1) ) + g2( in(n1+1:n2) );
            proxg = @(in,t) [ proxg1( in(1:n1), t );  proxg2(in(n1+1:n2)) ];
            proxgConj = @(in,s) proxConj( proxg, in, s );

            metric1 = @(in) 0.5 * norm( applyMFS( in ) - b )^2;
            metric2 = @metricSpiritRegViolation;
            metrics = { metric1, metric2 };
            metricNames = { 'sparsity', 'violation' };
            [img, objValues, mValues] = pdhgWLS( img0, proxf, proxgConj, 'A', applyA, 'f', f, 'g', g, ...
              'N', 1000, 'metrics', metrics, 'metricNames', metricNames, 'printEvery', 1, 'verbose', verbose );   %#ok<ASGLU>

          else

            % minimize || M F S x ||_2
            img0 = mri_reconRoemer( img0 );
            b = applyM( kData );
            [ img, lsqrFlag, lsqrRelRes, lsqrIter, lsqrResVec ] = lsqr( ...
              @apply_vMFS, b, tol, nMaxIter, [], [], img0(:) );   %#ok<ASGLU>
            img = reshape( img, sImg );
          end
  
        end
  
      else  % if numel( sMaps ) > 0
  
        if numel( support ) > 0
  
          if epsSupport > 0
            % minimize || M F x - b ||_2  subject to  || Pc x(c) ||_2^2 / ( Nc - 1 ) <= epsSupport for all c
  
            sampleMaskSlice = sampleMask(:,:,1);
            coilImgs = cell( 1, 1, nCoils );
            for coilIndx = 1 : nCoils
              coilImg0 = zeros( sImg );
              b = applyM( kData(:,:,coilIndx) );
              g = @(in) 0.5 * norm( applyF( in ) .* sampleMaskSlice - b, 'fro' )^2;
              FTMTb = applyF( kData, 'transp' );
              gGrad = @(in) applyF( applyF( in ) .* sampleMaskSlice, 'transp' ) - FTMTb;
              proxth = @(in) projectOntoBall( applyPC( in ), sqrt( epsSupport * (nSupportC - 1) ) );
              [coilImg,objectiveValues,relDiffs] = fista_wLS( coilImg0, g, gGrad, proxth );   %#ok<ASGLU>
              coilImgs{ 1, 1, coilIndx } = coilImg;
            end
            coilImgs = cell2mat( coilImgs );
            img = mri_reconRoemer( coilImgs );
  
          else
            % minimize || M F PcoilsT x - b ||_2
            b = applyM( kData );
            x0 = applyPcoils( coilImgs0 );
            [ coilImgsSupport, lsqrFlag, lsqrRelRes, lsqrIter, lsqrResVec ] = lsqr( ...
              @apply_vMFPcoilsT, b, tol, nMaxIter, [], [], x0(:) );   %#ok<ASGLU>
            coilImgs = applyPcoils( reshape( coilImgsSupport, nSupport, nCoils ), 'transp' );
            img = mri_reconRoemer( coilImgs );
          end
  
        else  % if numel( support ) > 0
          coilRecons = applyF( kData, 'inv' );
          img = mri_reconRoemer( coilRecons );
  
        end
  
      end
  
    else  % if nCoils > 1
  
      if numel( support ) > 0
  
        if epsSupport > 0
          % minimize || M F x - b ||_2  subject to  || Pc x ||_2^2 / ( Nc - 1 ) <= epsSupport
          % minimize (1/2) || M F x - b ||_2^2  subject to  || Pc x ||_2  <= sqrt( epsSupport * ( Nc - 1 ) )
          b = applyM( kData );
          g = @(in) 0.5 * norm( applyMF( in ) - b, 'fro' )^2;
          FTMTb = applyF( kData, 'transp' );
          %gGrad = @(in) applyMF( applyMF( in ), 'transp' ) - FTMTb;
          gGrad = @(in) applyF( applyF( in ) .* sampleMask, 'transp' ) - FTMTb;
          proxth = @(in) projectOntoBall( applyPC( in ), sqrt( epsSupport * (nSupportC - 1) ) );
          [img,objectiveValues,relDiffs] = fista_wLS( img0, g, gGrad, proxth );   %#ok<ASGLU>
  
        else
          if numel( optAlg ) == 0, optAlg = 'pocs'; end
          if strcmp( optAlg, 'pocs' )
            % minimize 1 subject to Pc x = 0  and  M F x = b
            proj1 = @(in) applyF( applyF( in ) .* ( 1 - sampleMask ) + kData, 'inv' );
            proj2 = @(in) in .* support;
            projections = { proj2, proj1 };
            img = pocs( img0, projections );

          elseif strcmp( optAlg, 'gd' )
            % minimize 0.5 * || M F PT x - b ||_2^2
            PFTMTb = applyMFPT( b, 'transp' );
            g = @(in) 0.5 * norm( applyMFPT(in) - b )^2;
            gGrad = @(in) applyMFPT( applyMFPT( in ), 'transp' ) - PFTMTb;
            [img,oVals,relDiffs] = gradDescent( applyP(img0), gGrad, 'g', g, 'useLineSearch', true );   %#ok<ASGLU>

          elseif strcmp( optAlg, 'lsqr' )
            % minimize || M F PT x - b ||_2
            % Note: DOESN'T WORK!!!  This problem is not appropriate for LSQR
            img0 = applyP( img0 );
            if nCoils == 1
              [ img, lsqrFlag, lsqrRelRes, lsqrIter, lsqrResVec ] = lsqr( ...
                @applyMFPT, b, tol, nMaxIter, [], [], img0(:) );   %#ok<ASGLU>
            else
              [ img, lsqrFlag, lsqrRelRes, lsqrIter, lsqrResVec ] = lsqr( ...
                @apply_vMFPT, b, tol, nMaxIter, [], [], img0(:) );   %#ok<ASGLU>
            end
            img = applyP( img, 'transp' );
          end
        end
  
      else
        img = applyF( kData, 'inv' );
      end
    end

  end

end


%% SUPPORT FUNCTIONS

function [w, epsSpiritReg] = findW( acr, wSize )
  %-- Find the interpolation coefficients w
  sACR = size( acr );
  nCoils = sACR( 3 );

  w = cell( 1, 1, 1, nCoils );
  epsSpiritReg = zeros( nCoils, 1 );
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
    epsSpiritReg( coilIndx ) = norm( A * wCoil - b )^2 / numel(b);
    wCoil = [ wCoil( 1 : pt2RemoveIndx - 1 ); 0; wCoil( pt2RemoveIndx : end ); ];
    w{1,1,1,coilIndx} = reshape( wCoil, [ wSize nCoils ] );
  end
  w = cell2mat( w );
  %gammas = mean( gammas ) / nCoils;
end

