
function [ img, objValues, mValues, residualErrs ] = mriRecon( kData, varargin )
  % img = mri_Recon( kData [, 'epsData', epsData, 'epsSupport', epsSupport, 'gamma', gamma, 
  %       'lambda', lambda, 'nOptIter', nOptIter, 'Psi', Psi 'sACR', sACR, 'sMaps', sMaps,
  %       'support', support, 'verbose', verbose, 'wSize', wSize ] )
  %
  % Inputs:
  % kData - an array of size M x N x C where C is the number of coils
  %         and uncollected data have value 0.
  %
  % Optional Inputs:
  % sMaps - a complex array of size M x N x C specifying the sensitivity of each coil
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
  p.addParameter( 'debug', false );
  p.addParameter( 'doChecks', false );
  p.addParameter( 'epsData', [], @isnonnegative );
  p.addParameter( 'epsSpiritReg', [], @isnonnegative );
  p.addParameter( 'epsSupport', 0, @isnonnegative );
  p.addParameter( 'lambda', [], @isnonnegative );
  p.addParameter( 'nOptIter', 2500, @ispositive );
  p.addParameter( 'noCS', false );
  p.addParameter( 'oPsi', true );
  p.addParameter( 'optAlg', [] );
  p.addParameter( 'Psi', [] );
  p.addParameter( 'sACR', [] );
  p.addParameter( 'sMaps', [], @isnumeric );
  p.addParameter( 'support', [] );
  p.addParameter( 'verbose', true );
  p.addParameter( 'wSize', [], @ispositive );
  p.parse( varargin{:} );
  debug = p.Results.debug;
  doChecks = p.Results.doChecks;
  epsData = p.Results.epsData;
  epsSpiritReg = p.Results.epsSpiritReg;
  epsSupport = p.Results.epsSupport;
  lambda = p.Results.lambda;
  nOptIter = p.Results.nOptIter;
  noCS = p.Results.noCS;
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

  if noCS == true
    cs = false;
  elseif ( numel( lambda ) > 0  &&  lambda > 0 )  ||  numel( epsData ) > 0
    if ( numel( lambda ) > 0  &&  lambda > 0 )  &&  numel( epsData ) > 0
      error( 'Cannot supply both lambda and epsData')
    end
    cs = true;
  else
    cs = false;
  end

  coilImgs0 = applyF( kData, 'inv' );
  if nCoils > 1
    img0 = mri_reconRoemer( coilImgs0, 'sMaps', sMaps );
  else
    img0 = coilImgs0;
  end

  spiritScaling = 1;
  if numel( wSize ) > 0
    if nCoils == 1, error( 'Cannot supply wSize with a single coil' ); end
    acr = cropData( kData, [ sACR(1) sACR(2) nCoils ] );
    if noCS == true
      spiritNormWeights = [];
      spiritNormWeightsACR = [];
    else
      spiritNormWeights = findSpiritNormWeights( kData );
      spiritNormWeightsACR = cropData( spiritNormWeights, [ sACR(1) sACR(2) nCoils ] );
    end
    if numel( epsSpiritReg ) == 0
      [w, epsSpiritReg] = findW( acr, wSize, spiritNormWeightsACR );
    elseif isscalar( epsSpiritReg )
      epsSpiritReg = epsSpiritReg * ones( nCoils, 1 );
      [w, ~] = findW( acr, wSize, spiritNormWeightsACR );
    else
      [w, ~] = findW( acr, wSize, spiritNormWeightsACR );
    end
    flipW = padData( flipDims( w, 'dims', [1 2] ), [ sImg nCoils nCoils ] );

    if numel( sMaps ) > 0
      normMFS = powerIteration( @applyMFS, rand( sImg ) );
      normSwWmIFS = powerIteration( @applySwWmIFS, rand( sImg ) );
      spiritScaling = normMFS / normSwWmIFS;
    end
  end

  if numel( Psi ) == 0
    wavSplit = makeWavSplit( sImg );
    Psi = @defaultPsi;
    if oPsi == false
      warning( 'The default value of oPsi being overwritten because default Psi is orthogonal' );
    end
    oPsi = true;
  end

  if numel( support ) > 0
    if any( sImg ~= size( support ) ), error( 'Support must be the size of the image' ); end
    nSupport = sum( support(:) );
    nSupportC = numel( support ) - nSupport;
  end

  if numel( sMaps ) > 0  &&  nCoils == 1
    error( 'Cannot supply sMaps with single coil data' );
  end

  b = applyM( kData );
  nb = numel( b );
  nbPerCoil = nb / nCoils;

  %-- Function definitions

  function out = applyA_SPIRiT( in, op )
    if nargin < 2 || strcmp( op, 'notransp' )
      out = A_SPIRiT( in );
    else
      out = AH_SPIRiT( reshape( in, size( kData ) ) );
    end
    out = out(:);
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

  function out = applyMc( in, op )
    if nargin < 2  ||  strcmp( op, 'notransp' )
      out = in( sampleMask == 0 );
    else
      out = zeros( sKData );
      out( sampleMask == 0 ) = in;
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

  function out = applyMFSPre( in, op )
    if nargin < 2  ||  strcmp( op, 'notransp' )
      out = applyMFS( applyPre( in ) );
    else
      out = applyPre( applyMFS( in, op ), op );
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
      for cI = 1 : nCoils
        out(:,cI) = applyP( in(:,:,cI) );
      end
    else
      out = zeros( sKData );
      for cI = 1 : nCoils
        out(:,:,cI) = applyP( in(:,cI), op );
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
      for cI = 1 : nCoils
        out(:,cI) = applyPC( in(:,:,cI) );
      end
    else
      out = zeros( sKData );
      for cI = 1 : nCoils
        out(:,:,cI) = applyPC( in(:,cI), op );
      end
    end
  end

  S_norms = LpNorms( sMaps, 2, 3 );
  function out = applyPre( in, op )
    if nargin < 2  ||  strcmp( op, 'notransp' )
      out = in ./ S_norms;
    elseif strcmp( op, 'transp' )
      out = in ./ conj( S_norms );
    elseif strcmp( op, 'inv' )
      out = in .* S_norms;
    elseif strcmp( op, 'invTransp' )
      out = in .* conj( S_norms );
    end
  end

  function out = applyS( in, op )
    if nargin < 2 || strcmp( op, 'notransp' )
      out = bsxfun( @times, sMaps, in );
    else
      out = sum( conj(sMaps) .* in, 3 );
    end
  end

  function out = applySw( in, op )   %#ok<INUSD>
    % Apply the spirit norm weights
    % Note that since the weights are real, Sw is self adjoint
    out = zeros( size( in ) );
    for cI = 1 : nCoils
      %out(:,:,cI) = 1d5 * spiritNormWeights(:,:,cI) .* in(:,:,cI);
      out(:,:,cI) = spiritScaling * spiritNormWeights(:,:,cI) .* in(:,:,cI);
    end
  end

  function out = applySwWmI( in, op )
    % apply Sw .* ( W minus I )
    if nargin < 2 || strcmp( op, 'notransp' )
      out = applySw( applyWmI( in ) );
    else
      SwHIn = applySw( in, 'transp' );
      out = applyWmI( SwHIn, 'transp' );
    end
  end

  function out = applySwWmIFS( in, op )
    % apply Sw .* ( W minus I )
    if nargin < 2 || strcmp( op, 'notransp' )
      out = applySw( applyWmI( applyFS( in ) ) );
    else
      SwHIn = applySw( in, 'transp' );
      out = applyFS( applyWmI( SwHIn, 'transp' ), 'transp' );
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
      out = applyWmI( FSin );
    else
      out = applyFS( applyWmI( in, 'transp' ), 'transp' );
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

  bPerCoil = reshape( b, [], nCoils );
  function out = indicatorEpsData( in )
    % indicator( || x(c) - b(c) ||_2  <= ballRadius(c) )
    out = 0;    
    in = reshape( in, [], nCoils );
    for cI = 1 : nCoils
      out = indicatorFunction( norm( in - bPerCoil(:,cI) )^2 / nbPerCoil, [ 0 epsData(cI) ] );
      if out ~= 0, break; end
    end
  end

  epsBallRadiuses = sqrt( epsData * nbPerCoil );
  function out = projectOntoEpsBalls( in )
    % proximal operator of indicator( || x(c) - b(c) ||_2  <= ballRadius(c) )
    % In can either be a vector of the size of b, or it can be a vector of the size of kData
    % If the size of kData, the projection is only done for those values that were sampled

    allOut = [];
    if numel( in ) == numel( kData )
      allOut = in;
      in = in( sampleMask == 1 );
    end

    sOut = size( in );
    in = reshape( in, [], nCoils );
    out = zeros( size( in, 1 ), nCoils );
    for cI = 1 : nCoils
      out(:,cI) = projectOntoBall( in(:,cI) - bPerCoil(:,cI), epsBallRadiuses(cI) ) + bPerCoil(:,cI);
    end
    out = reshape( out, sOut );

    if numel( allOut ) > 0
      allOut( sampleMask == 1 ) = out;
      out = allOut;
    end
  end

  function out = metricEpsDataViolation( in )
    MFSin = applyMFS( in );
    MFSin = reshape( MFSin, [], nCoils );

    violations = zeros( nCoils, 1 );
    for cI = 1 : nCoils
      cProj = projectOntoBall( MFSin(:,cI) - bPerCoil(:,cI), epsBallRadiuses(cI) ) + bPerCoil(:,cI);
      violations(cI) = norm( cProj - MFSin(:,cI) );
    end
    out = max( violations );
  end

  if numel( support ) > 0  &&  numel( epsSupport ) > 0
    outsideSupportBallRadius = sqrt( epsSupport * nSupportC );
  end

  function out = indicatorOutsideSupport( in )
    out = indicatorFunction( norm( applyPC( in ) ), [ 0, outsideSupportBallRadius] );
  end

  function out = projOutsideSupportOntoBall( in )
    out = in;
    out( support == 0 ) = projectOntoBall( applyPC( in ), outsideSupportBallRadius );
  end

  function out = metricOutsideSupportViolation( in )
    out = max( norm( applyPC( in ) ) - outsideSupportBallRadius, 0 );
  end

  if numel( wSize ) > 0
    ballRadiusesSpiritReg = sqrt( epsSpiritReg * nImg );
  end

  function out = indicatorSpiritReg( x )
    % indicator( || x(:,:,coilIndx) ||_2^2 / nImg <= epsSpiritReg(coilIndx) ) is equivalent to
    % indicator( || x(:,:,coilIndx) ||_Fro <= sqrt( epsSpiritReg( coilIndx ) * nImg ) )
    xNorms = LpNorms( reshape( x, [], nCoils ), 2, 1 );
    out = 0;
    if any( xNorms > ballRadiusesSpiritReg )
      out = Inf;
    end
  end

  function out = proxIndSpiritReg( in, t )   %#ok<INUSD>
    % indicator( || in(:,:,coilIndx) ||_2^2 / nImg <= epsSpiritReg(coilIndx) ) is equivalent to
    % indicator( || in(:,:,coilIndx) ||_Fro <= sqrt( epsSpiritReg( coilIndx ) * nImg ) )
    in = reshape( in, sKData );
    out = zeros( sKData );
    for cI = 1 : nCoils
      out(:,:,cI) = projectOntoBall( in(:,:,cI), ballRadiusesSpiritReg(cI), 'fro' );
    end
  end

  function out = metricSpiritRegViolation( in )
    WmIFSin = applyWmIFS( in );
    violations = zeros( nCoils, 1 );
    for cI = 1 : nCoils
      violations(cI) = max( norm( WmIFSin(:,:,cI), 'fro' ) - ballRadiusesSpiritReg(cI), 0 );
    end
    out = max( violations );
  end

  function out = metricWeightedSpiritRegViolation( in )
    SwWmIFSin = applySw( applyWmIFS( in ) );
    violations = zeros( nCoils, 1 );
    for cI = 1 : nCoils
      violations(cI) = max( norm( SwWmIFSin(:,:,cI), 'fro' ) - ballRadiusesSpiritReg(cI), 0 );
    end
    out = max( violations );
  end


  %-- Concatenated functions

  function out = concat_MFS_I( in, op )
    if nargin < 2  ||  strcmp( op, 'notransp' )
      out = [ reshape( applyMFS(in), [], 1 );  in(:) ];
    else
      out = applyMFS( in(1:nb), op ) + reshape( in(nb+1:end), sImg );
    end
  end

  function out = concat_MFS_Pc( in, op )
    if nargin < 2  ||  strcmp( op, 'notransp' )
      out = [ reshape( applyMFS(in), [], 1 );  applyPC( in ); ];
    else
      out = applyMFS( in(1:nb), op ) + reshape( applyPC( in(nb+1:end), 'transp' ), sImg );
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

  function out = concat_MFS_SwWmIFS( in, op )
    if nargin < 2 || strcmp( op, 'notransp' )
      FSin = applyFS( in );
      MFSin = applyM( FSin );
      SwWmIFSin = applySwWmI( FSin );
      out = [ MFSin(:); SwWmIFSin(:); ];
    else
      n1 = nb;
      n2 = n1 + nKData;
      MTin1 = reshape( applyM( in( 1 : n1 ), op ), sKData );
      SwHWHmIHin2 = reshape( applySwWmI( reshape( in( n1+1 : n2 ), sKData ), op ), sKData );
      out = applyFS( MTin1 + SwHWHmIHin2, op );
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

  function out = concat_MFS_SwWmIFS_I( in, op )
    if nargin < 2 || strcmp( op, 'notransp' )
      FSin = applyFS( in );
      MFSin = applyM( FSin );
      SwWmIFSin = applySw( applyWmI( FSin ) );
      out = [ MFSin(:); SwWmIFSin(:); in(:); ];
    else
      n1 = nb;
      n2 = n1 + nKData;
      n3 = n2 + nImg;
      MTin1 = reshape( applyM( in( 1 : n1 ), op ), sKData );
      WHmIHSwTin2 = reshape( applyWmI( applySw( reshape( in( n1+1 : n2 ), sKData ), op ), op ), sKData );
      in3 = reshape( in( n2+1 : n3 ), sImg );
      out = applyFS( MTin1 + WHmIHSwTin2, op ) + in3;
    end
  end

  function out = concat_Psi_MFS( in, op )
    if nargin < 2  ||  strcmp( op, 'notransp' )
      out = [ reshape( Psi(in), [], 1 ); ...
              reshape( applyMFS(in), [], 1 ); ];
    else
      n1 = nPsi;
      n2 = nPsi + nb;
      out = Psi( in(1:nPsi), op ) + applyMFS( in(n1:n2), op );
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

  function out = concat_Psi_MFS_SwWmIFS( in, op )
    if nargin < 2 || strcmp( op, 'notransp' )
      PsiIn = Psi( in );
      FSin = applyFS( in );
      MFSin = applyM( FSin );
      SwWmIFSin = applySwWmI( FSin );
      out = [ PsiIn(:); MFSin(:); SwWmIFSin(:); ];
    else
      n0 = nPsi;
      n1 = n0 + nb;
      n2 = n1 + nKData;
      PsiIn0 = Psi( reshape( in( 1 : n0 ), sKData), op );
      MTin1 = reshape( applyM( in( n0 : n1 ), op ), sKData );
      WHmIHSwTin3 = reshape( applySwWmI( reshape( in( n1+1 : n2 ), sKData ), op ), sKData );
      out = applyFS( PsiIn0 + MTin1 + WHmIHSwTin3, op );
    end
  end


  %-- Vectorized functions

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

  function out = apply_vMFSPre( in, op )
    if nargin < 2  ||  strcmp( op, 'notransp' )
      in = reshape( in, sImg );
      out = applyMFSPre( in );
    else
      out = applyMFSPre( in, op );
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

  function out = v_proxIndSpiritReg( in, t )
    out = proxIndSpiritReg( in, t );
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

    [chk_MF,errChk_MF] = checkAdjoint( rand( sKData ) + 1i * rand( sKData ), @applyMF );
    if chk_MF == true
      disp( 'Check of applyMF passed' );
    else
      error([ 'Check of applyMF Adjoint failed with error ', num2str(errChk_MF) ]);
    end

    [chk_MFS,errChk_MFS] = checkAdjoint( rand( sKData ) + 1i * rand( sKData ), @applyMFS );
    if chk_MFS == true
      disp( 'Check of applyMFS passed' );
    else
      error([ 'Check of applyMFS Adjoint failed with error ', num2str(errChk_MFS) ]);
    end

    [chk_Pre,errChk_Pre] = checkAdjoint( rand( sKData ) + 1i * rand( sKData ), @applyPre );
    if chk_Pre == true
      disp( 'Check of applyPre passed' );
    else
      error([ 'Check of applyPre Adjoint failed with error ', num2str(errChk_Pre) ]);
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
      [chkSw,errChkSw] = checkAdjoint( rand( sKData ) + 1i * rand( sKData ), @applySw );
      if chkSw == true
        disp( 'Check of Sw Adjoint passed' );
      else
        error([ 'Check of Sw Adjoint failed with error ', num2str(errChkSw) ]);
      end

      [chkW,errChkW] = checkAdjoint( rand( sKData ) + 1i * rand( sKData ), @applyW );
      if chkW == true
        disp( 'Check of W minus I Adjoint passed' );
      else
        error([ 'Check of W minus I Adjoint failed with error ', num2str(errChkW) ]);
      end

      [chkWmI,errChkWmI] = checkAdjoint( rand( sKData ) + 1i * rand( sKData ), @applyWmI );
      if chkWmI == true
        disp( 'Check of W minus I Adjoint passed' );
      else
        error([ 'Check of W minus I Adjoint failed with error ', num2str(errChkWmI) ]);
      end

      [chkConcat_MFS_SwWmIFS,errChkConcat_MFS_SwWmIFS] = checkAdjoint( rand( sKData(1:2) ) + 1i * rand( sKData(1:2) ), @concat_MFS_SwWmIFS );
      if chkConcat_MFS_SwWmIFS == true
        disp( 'Check of concat_MFS_SwWmIFS Adjoint passed' );
      else
        error([ 'Check of concat_MFS_SwWmIFS Adjoint failed with error ', num2str(errChkConcat_MFS_SwWmIFS) ]);
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
        [chkWmIFS,errChkWmIFS] = checkAdjoint( rand( sImg ) + 1i * rand( sImg ), @applyWmIFS );
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

        [ chk_MFS_I, errChk_MFS_I ] = checkAdjoint( rand( sImg ) + 1i * rand( sImg ), @concat_MFS_I );
        if chk_MFS_I == true
          disp( 'Check of concat_MFS_I adjoint passed' );
        else
          error([ 'Check of concat_MFS_I adjoint failed with error ', num2str(errChk_MFS_I) ]);
        end
      end
    end
  end  % if doChecks


  %-- Reconstruct the image

  tol = 1d-6;

  if cs == true

    nPsiCoils = numel( Psi( coilImgs0 ) );
    nPsi = numel( Psi( img0 ) );

    if nCoils > 1
  
      if numel( sMaps ) > 0

        if numel( support ) > 0

          if epsSupport > 0

            if numel( wSize ) > 0
              % minimize || Psi x ||_1
              % subject to  || M F S(c) x - b(c) ||_2^2 / nbPerCoil <= epsData(c)
              %        and  || Sw(c) ((W-I) F S x)(:,:,c) ||_Fro <= sqrt( epsSpiritReg(c) * nImg ) for all c
              %        and  || Pc x ||_2^2 / ( nImg - Ns - 1 ) <= epsSupport

              g1 = @( in ) indicatorEpsData( in );
              proxg1 = @(in,t) projectOntoEpsBalls( in );

              if oPsi == true
                applyA = @(in,op) concat_MFS_SwWmIFS_I( in, op );
                f = @(in) norm( reshape( Psi( in ), [], 1 ), 1 );
                proxf = @(in,t) proxCompositionAffine( @proxL1Complex, in, Psi, 0, 1, t );
                g2 = @indicatorSpiritReg;
                proxg2 = @v_proxIndSpiritReg;
                g3 = @(in) indicatorFunction( ...
                  norm( in, 2 )^2 / ( nImg - nSupport - 1 ), [ 0 epsSupport ] );
                proxg3 = @(in,t) projOutsideSupportOntoBall( in );
                n1 = nb;
                n2 = n1 + nKData;
                n3 = n2 + nImg;
                g = @(in) g1( in(1:n1) ) + g2( in(n1+1:n2) ) + g3( in(n2+1:n3) );
                proxg = @(in,t) [ proxg1( in(1:n1), t );  ...
                                  proxg2( in(n1+1:n2), t ); ...
                                  proxg3( in(n2+1:n3) ); ];
                proxgConj = @(in,s) proxConj( proxg, in, s );
    
                PsiImg0 = Psi( img0 );
                tau0 = mean( abs( PsiImg0(:) ) );

                metric1 = @(in) 0.5 * norm( applyMFS( in ) - b )^2;
                metric2 = @metricSpiritRegViolation;
                metrics = { metric1, metric2 };
                metricNames = { 'sparsity', 'spirit violation' };
                [img, objValues, mValues] = pdhgWLS( img0, proxf, proxgConj, 'A', applyA, 'f', f, 'g', g, ...
                  'N', nOptIter, 'tau', tau0, 'metrics', metrics, 'metricNames', metricNames, ...
                  'printEvery', 10, 'verbose', verbose );

              else
                applyA = @(in,op) concat_Psi_MFS_SwWmIFS( in, op );
                f = @(in) indicatorFunction( ...
                  norm( in, 2 )^2 / ( nImg - nSupport - 1 ), [ 0 epsSupport ] );
                proxf = @(in,t) projOutsideSupportOntoBall( in );
                g0 = @(in) norm( in(:), 1 );
                proxg0 = @proxL1Complex;
                g2 = @indicatorSpiritReg;
                proxg2 = @v_proxIndSpiritReg;
                n0 = nKData;
                n1 = n0 + nb;
                n2 = n1 + nKData;
                g = @(in) g0( in(1:n0) ) + g1( in(n0:n1) ) + g2( in(n1+1:n2), t );
                proxg = @(in,t) [ proxg0( in(1:n0), t ); ...
                                  proxg1( in(n0+1:n1), t ); ...
                                  proxg2( in(n1+1:n2), t ) ];
                proxgConj = @(in,s) proxConj( proxg, in, s );

                metric1 = @(in) 0.5 * norm( applyMFS( in ) - b )^2;
                metric2 = @metricSpiritRegViolation;
                metrics = { metric1, metric2 };
                metricNames = { 'sparsity', 'spirit violation' };
                [img, objValues, mValues] = pdhgWLS( img0, proxf, proxgConj, 'A', applyA, 'f', f, 'g', g, ...
                  'N', nOptIter, 'metrics', metrics, 'metricNames', metricNames, 'printEvery', 10, 'verbose', verbose );   %#ok<ASGLU>
              end
            else

              if numel( epsData ) > 0
                % minimize || Psi x ||_1
                % subject to || M F S(c) x - b(c) ||_2^2 / nbPerCoil <= epsData(c)
                % and  || Pc x ||_2^2 / ( nImg - Ns - 1 ) <= epsSupport
  
                g1 = @( in ) indicatorEpsData( in );
                proxg1 = @(in,t) projectOntoEpsBalls( in );
                if oPsi == true
                  f = @(in) norm( reshape( Psi( in ), [], 1 ), 1 );
                  proxf = @(in,t) proxCompositionAffine( @proxL1Complex, in, Psi, 0, 1, t );
                  applyA = @concat_MFS_I;
                  g2 = @(in) indicatorFunction( ...
                    norm( in( support == 0 ), 2 )^2 / ( nImg - nSupport ), [ 0 epsSupport ] );
                  g = @(in) g1( in(1:nb) ) + g2( in(nb+1:end) );
                  proxg2 = @(in,t) projOutsideSupportOntoBall( in );
                  proxg = @(in,t) [ proxg1( in(1:nb), t );  proxg2( in(nb+1:end), t ); ];
                  proxgConj = @(in,s) proxConj( proxg, in, s );

                  PsiImg0 = Psi( img0 );
                  tau0 = mean( abs( PsiImg0(:) ) );

                  [img, objValues, mValues] = pdhgWLS( img0, proxf, proxgConj, 'A', applyA, ...
                    'f', f, 'g', g, 'tau', tau0, 'printEvery', 10, 'verbose', verbose );
  
                else  % if oPsi == true
                  applyA = @concat_Psi_MFS;
                  f = @(in) indicatorFunction( ...
                    norm( in( support == 0 ), 2 )^2 / ( nImg - nSupport ), [ 0 epsSupport ] );
                  proxf = @(in,t) projOutsideSupportOntoBall( in );
                  n1 = nImg;
                  n2 = n1 + nb;
                  g0 = @(in) norm( reshape( in, [], 1 ), 1 );
                  g = @(in) g0( in(1:n1) ) + g1( in(n1+1:n2) );
                  proxg0 = @proxL1Complex;
                  proxg = @(in,t) [ proxg0( in(1:n1), t ); ...
                                    proxg1( in(n1+1:n2), t ); ];
                  proxgConj = @(in,s) proxConj( proxg, in, s );
                  [img, objValues, mValues] = pdhgWLS( img0, proxf, proxgConj, 'A', applyA, ...
                    'f', f, 'g', g, 'printEvery', 10, 'verbose', verbose );
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
              % subject to  || M F S(c) x - b ||_2^2 / nbPerCoil <= epsData(c)
              %        and  || Sw ((W-I) F S x)(:,:,c) ||_Fro <= sqrt( epsSpiritReg(c) * nImg ) for all c

              g1 = @( in ) indicatorEpsData( in );
              proxg1 = @(in,t) projectOntoEpsBalls( in );
              g2 = @indicatorSpiritReg;
              proxg2 = @v_proxIndSpiritReg;

              metric1 = @(in) norm( in(:), 1 );
              metric2 = @metricEpsDataViolation;
              metric3 = @metricSpiritRegViolation;
              metrics = { metric1, metric2, metric3 };
              metricNames = { 'sparsity', 'dcViolation', 'spiritViolation' };

              if oPsi == true
                applyA = @(in,op) concat_MFS_SwWmIFS( in, op );
                f = @(in) norm( reshape( Psi( in ), [], 1 ), 1 );
                proxf = @(in,t) proxCompositionAffine( @proxL1Complex, in, Psi, 0, 1, t );
                % || M F S x - b ||_2^2 / (nb-1) <= epsData  <==>  || M F S x - b ||_2 <= sqrt( epsData * (nb-1) )
                n1 = nb;
                n2 = n1 + nKData;
                g = @(in) g1( in(1:n1) ) + g2( in(n1+1:n2) );
                proxg = @(in,t) [ proxg1( in(1:n1), t );  proxg2(in(n1+1:n2), t) ];
                proxgConj = @(in,s) proxConj( proxg, in, s );
    
                PsiImg0 = Psi( img0 );
                tau0 = mean( abs( PsiImg0(:) ) );

                [img, objValues, mValues] = pdhgWLS( img0, proxf, proxgConj, 'A', applyA, 'f', f, 'g', g, ...
                  'N', nOptIter, 'tau', tau0, 'metrics', metrics, 'metricNames', metricNames, ...
                  'printEvery', 10, 'verbose', verbose );

              else
                applyA = @(in,op) concat_Psi_MFS_SwWmIFS( in, op );
                f = [];
                proxf = [];
                g0 = @(in) norm( in, 1 );
                proxg0 = @proxL1Complex;
                n0 = nKData;
                n1 = n0 + nb;
                n2 = n1 + nKData;
                g = @(in) g0( in(1:n0) ) + g1( in(n0+1:n1) ) + g2( in(n1+1:n2) );
                proxg = @(in,t) [ proxg0( in(n0:n1), t ); ...
                                  proxg1( in(n0+1:n1), t ); ...
                                  proxg2( in(n1+1:n2), t) ];
                proxgConj = @(in,s) proxConj( proxg, in, s );

                [img, objValues, mValues] = pdhgWLS( img0, proxf, proxgConj, 'A', applyA, 'f', f, 'g', g, ...
                  'N', nOptIter, 'metrics', metrics, 'metricNames', metricNames, 'printEvery', 10, 'verbose', verbose );   %#ok<ASGLU>
              end
            end

          else % if numel( wSize ) > 0

            % Parallel imaging with compressed sensing (PICS)

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
                if debug == true
                  [ img, objValues, relDiffs] = fista_wLS( img0, f, fGrad, proxg, 'h', g, 't0', t0, ...
                    'printEvery', 50, 'verbose', verbose );   %#ok<ASGLU>
                else
                  img = fista_wLS( img0, f, fGrad, proxg, 'h', g );
                end
  
              else  % if numel( lambda ) > 0  &&  lambda > 0
                % minimize || Psi x ||_1  subject to  || M F S(c) x - b ||_2^2 / ( nbPerCoil - 1 ) <= epsData(c)

                applyA = @applyMFS;
                f = @(in) norm( vPsi( in ), 1 );
                proxf = @(in,t) proxCompositionAffine( @proxL1Complex, in, Psi, [], 1, t );
                g = @( in ) indicatorEpsData( in );
                proxg = @(in,t) projectOntoEpsBalls( in );
                proxgConj = @(in,s) proxConj( proxg, in, s );
  
                metric1 = @(in) 0.5 * norm( applyMFS( in ) - b )^2;
                metric2 = @metricEpsDataViolation;
                metrics = { metric1, metric2 };
                metricNames = { 'data consistency', 'violation' };

                PsiImg0 = Psi( img0 );
                tau0 = mean( abs( PsiImg0(:) ) );

                [img, objValues, mValues] = pdhgWLS( img0, proxf, proxgConj, 'A', applyA, 'f', f, 'g', g, ...
                  'tau', tau0, 'metrics', metrics, 'metricNames', metricNames, 'printEvery', 10, 'verbose', verbose );   %#ok<ASGLU>
  
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
            % minimize || Psi x ||_{2,1} subject to || M F x - b ||_2^2 / nbPerCoil < epsData(c)
            % and || Pc x ||_2^2 / nOutsideSupport < epsSupport
            error( 'Not yet implemented' );

          else
            error( 'Also not yet implemented' );

          end

        else % if numel( support ) > 0

          if numel( epsData ) > 0
            % minimize || Psi x ||_{2,1} subject to || M F x - b ||_2^2 / nbPerCoil < epsData(c)

            if oPsi == true
              applyA = @applyMF;
              f = @(in) normL2L1( Psi( in ) );
              proxf = @(in,t) proxCompositionAffine( @proxL2L1, in, Psi, 0, 1, t );
              g = @( in ) indicatorEpsData( in );
              proxg = @(in,t) projectOntoEpsBalls( in );
              proxgConj = @(in,s) proxConj( proxg, in, s );
              [y, objValues, mValues] = pdhgWLS( img0, proxf, proxgConj, 'A', applyA, ...
                'beta', 1, 'f', f, 'g', g, 'N', N, 'dsc', true, ...
                'printEvery', 10, 'metrics', metrics, 'metricNames', metricNames, 'verbose', verbose );   %#ok<ASGLU>

            else  % if oPsi == true
              applyA = @concat_Psi_MF;
              f = [];
              proxf = [];
              g1 = @(in) normL2L1( in );
              g2 = @( in ) indicatorEpsData( in );
              g = @(in) g1( reshape(in(1:nKData), sKData) ) + g2( in(nKData+1:end) );
              proxg1 = @(in,t) proxL2L1( reshape( in, sKData ), t );
              proxg2 = @(in,t) projectOntoEpsBalls( in );
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
          error( 'Not yet implemented' );

        else  % if epsSupport > 0
          error( 'Also not yet implemented' );

        end

      else  % if numel( support ) > 0

        if numel( lambda ) > 0  &&  lambda > 0
          % minimize (1/2) || M F x - b ||_2^2 + lambda || Psi x ||_1

          if oPsi == true
            g = @(in) 0.5 * norm( applyMF( in ) - b, 2 )^2;
            FTMTb = applyF( kData, 'transp' );
            gGrad = @(in) applyF( applyF( in ) .* sampleMask, 'transp' ) - FTMTb;
            proxth = @(in,t) proxCompositionAffine( @proxL1Complex, in, Psi, [], 1, t * lambda );
            if debug == true
              img = fista_wLS( img0, g, gGrad, proxth, 'printEvery', 50, 'verbose', true );
            else
              img = fista_wLS( img0, g, gGrad, proxth );
            end

          else
          end

        else
          % minimize || Psi x ||_1  subject to  || M F x - b ||_2^2 / nbPerCoil <= epsData(c)
          % if A = M F the A AT = M F FT MT = a scaled identity matrix

          if oPsi == true
            % Solve this problem with ADMM
            error( 'Not yet implemented' );
          else
            error( 'Also not yet implemented' );
          end
        end

      end

    end

  else  % if cs == true

    if nCoils > 1

      if numel( sMaps ) > 0

        if numel( support ) > 0

          if epsSupport > 0

            if numel( wSize ) > 0  % do SPIRiT Regularization
              % minimize || M F S x - b ||_2 
              % subject to || Pc x ||_2^2 / ( Nc - 1 ) <= epsSupport
              % and || Sw ((W-I) F S x)(:,:,c) ||_Fro <= sqrt( epsSpiritReg(c) * nImg ) for all c

              applyA = @(in,op) concat_MFS_SwWmIFS( in, op );
              n1 = nb;
              n2 = n1 + nKData;
              f = @indicatorOutsideSupport;
              proxf = @(in,t) projOutsideSupportOntoBall( in );
              g1 = @( x ) 0.5 * norm( x - b )^2;
              proxg1 = @(in,t) proxL2Sq( in, t, b );
              g2 = @indicatorSpiritReg;
              proxg2 = @v_proxIndSpiritReg;
              g = @(in) g1( in(1:n1) ) + g2( in(n1+1:n2) );
              proxg = @(in,t) [ proxg1( in(1:n1), t );  proxg2(in(n1+1:n2), t) ];
              proxgConj = @(in,s) proxConj( proxg, in, s );

              metric1 = @(in) 0.5 * norm( applyMFS( in ) - b )^2;
              metric2 = @metricSpiritRegViolation;
              metric3 = @metricOutsideSupportViolation;
              metrics = { metric1, metric2, metric3 };
              metricNames = { 'sparsity', 'spirit violation', 'support violation' };
              if debug == true
                [img, objValues, mValues] = pdhgWLS( img0, proxf, proxgConj, 'A', applyA, 'f', f, 'g', g, ...
                  'N', nOptIter, 'metrics', metrics, 'metricNames', metricNames, 'verbose', verbose );   %#ok<ASGLU>
              else
                img = pdhgWLS( img0, proxf, proxgConj, 'A', applyA, 'f', f, 'g', g, ...
                  'N', nOptIter, 'metrics', metrics, 'metricNames', metricNames, 'printEvery', 10, 'verbose', verbose );
              end

            else
              % minimize || M F S x - b ||_2 subject to || Pc x ||_2^2 / ( Nc - 1 ) <= epsSupport
              g = @(in) 0.5 * norm( applyMFS( in ) - b, 2 )^2;
              SHFHMTb = applyFS( kData, 'transp' );
              %gGrad = @(in) applyMFS( applyMFS( in ), 'transp' ) - SHFHMTb;
              gGrad = @(in) applyFS( applyFS( in ) .* sampleMask, 'transp' ) - SHFHMTb;

              proxth = @(in,t) projOutsideSupportOntoBall( in );

              if numel( optAlg ) == 0, optAlg = 'projGrad'; end

              if debug == true
                if strcmp( optAlg, 'projGrad' )
                  [img,objValues] = projSubgrad( img0, gGrad, proxth, 'g', g, 'N', nOptIter, 'verbose', verbose );   %#ok<ASGLU>
                elseif strcmp( optAlg, 'fista_wLS' )
                  h = @(in) indicatorOutsideSupport( in );
                  [img,ojbValues,relDiffs] = fista_wLS( img0, g, gGrad, proxth, 'N', 200, 'h', h, ...
                    'printEvery', 10, 'verbose', verbose );   %#ok<ASGLU>
                else
                  error( 'Unrecognized optimization algorithm' );
                end

              else  % if debug == true
                if strcmp( optAlg, 'projGrad' )
                  img = projSubgrad( img0, gGrad, proxth, 'g', g );
                elseif strcmp( optAlg, 'fista_wLS' )
                  img = fista_wLS( img0, g, gGrad, proxth, 'N', nOptIter, 'verbose', verbose );
                else
                  error( 'Unrecognized optimization algorithm' );
                end
              end
              
            end
  
          else  % if epsSupport > 0

            % minimize || M F S PT x - b ||_2
            img0 = applyP( img0 );
            img0(:) = 0;
            b = applyM( kData );
            [ img, lsqrFlag, lsqrRelRes, lsqrIter, lsqrResVec ] = lsqr( ...
              @apply_vMFSPT, b, tol, nOptIter, [], [], img0(:) );   %#ok<ASGLU>
            img = applyP( img, 'transp' );

          end

        else  % if numel( support ) > 0

          if numel( wSize ) > 0  % do SPIRiT Regularization
            % minimize (1/2) || M F S x - b ||_2^2  
            % subject to || Sw(c) ((W-I) F S x)(:,:,c) ||_Fro^2 / nImg <= epsSpiritReg(c) for all c

            applyA = @(in,op) concat_MFS_SwWmIFS( in, op );
            n1 = nb;
            n2 = n1 + nKData;
            f = [];
            proxf = [];
            g1 = @( x ) 0.5 * norm( x - b )^2;
            proxg1 = @(in,t) proxL2Sq( in, t, b );
            g2 = @indicatorSpiritReg;
            proxg2 = @v_proxIndSpiritReg;
            g = @(in) g1( in(1:n1) ) + g2( in(n1+1:n2) );
            proxg = @(in,t) [ proxg1( in(1:n1), t );  proxg2(in(n1+1:n2), t) ];
            proxgConj = @(in,s) proxConj( proxg, in, s );

            metricDC = @( in ) 0.5 * norm( applyMFS( in ) - b )^2;
            metrics = { metricDC, @metricWeightedSpiritRegViolation };
            metricNames = { 'data consistency', 'violation' };
            [img, objValues, mValues] = pdhgWLS( img0, proxf, proxgConj, 'A', applyA, 'f', f, 'g', g, ...
              'N', nOptIter, 'metrics', metrics, 'metricNames', metricNames, 'printEvery', 10, 'verbose', verbose );

          else  % if numel( wSize ) > 0
            % minimize || M F S x - b ||_2

            if numel( optAlg ) == 0, optAlg = 'lsqr'; end

            if strcmp( optAlg, 'gd' )
              g = @(in) 0.5 * norm( applyMFS(in) - b )^2;
              %SHFHMTb = applyMFS( b, 'transp' );
              SHFHMTb = applyFS( kData, 'transp' );
              %gGrad = @(in) applyMFS( applyMFS( in ), 'transp' ) - SHFHMTb;
              gGrad = @(in) applyFS( applyFS( in ) .* sampleMask, 'transp' ) - SHFHMTb;              
              if debug
                [img,oVals,relDiffs] = gradDescent( img0, gGrad, 'g', g, 'useLineSearch', true, 'N', nOptIter, ...
                  'verbose', true, 'printEvery', 20 );   %#ok<ASGLU>
              else
                [img,oVals,relDiffs] = gradDescent( img0, gGrad, 'g', g, 'useLineSearch', true, 'N', nOptIter );   %#ok<ASGLU>
              end

            elseif strcmp( optAlg, 'lsqr' )

              usePreconditioner = false;
              if usePreconditioner == true
                [ y, lsqrFlag, lsqrRelRes, lsqrIter, lsqrResVec ] = lsqr( ...
                  @apply_vMFSPre, b, tol, nOptIter, [], [], img0(:) );   %#ok<ASGLU>
                img = applyPre( reshape( y, sImg ), 'inv' );
              else
                [ img, lsqrFlag, lsqrRelRes, lsqrIter, lsqrResVec ] = lsqr( ...
                  @apply_vMFS, b, tol, nOptIter, [], [], img0(:) );   %#ok<ASGLU>
                img = reshape( img, sImg );
              end
            end

            residualErrs = zeros( nCoils, 1 );
            for coilIndx = 1 : nCoils
              tmp = applyF( sMaps(:,:,coilIndx) * img );
              residualErrs( coilIndx ) = norm( tmp( sampleMask(:,:,coilIndx) == 1 ) - bPerCoil(:,coilIndx) )^2 / ...
                numel( bPerCoil(:,coilIndx) );
            end
            objValues = [];
            mValues = [];
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
              coilImg = fista_wLS( coilImg0, g, gGrad, proxth );
              coilImgs{ 1, 1, coilIndx } = coilImg;
            end
            coilImgs = cell2mat( coilImgs );
            img = mri_reconRoemer( coilImgs );

          else
            % minimize || M F PcoilsT x - b ||_2
            b = applyM( kData );
            x0 = applyPcoils( coilImgs0 );
            [ coilImgsSupport, lsqrFlag, lsqrRelRes, lsqrIter, lsqrResVec ] = lsqr( ...
              @apply_vMFPcoilsT, b, tol, nOptIter, [], [], x0(:) );   %#ok<ASGLU>
            coilImgs = applyPcoils( reshape( coilImgsSupport, nSupport, nCoils ), 'transp' );
            img = mri_reconRoemer( coilImgs );
          end

        else  % if numel( support ) > 0

          if numel( wSize ) > 0

            if numel( epsData ) > 0  &&  max( abs( epsData ) ) ~= 0
              % SPIRiT with soft consistency
              % minimize (1/2) || ( W - I ) k ||_2^2
              % subject to || M k(:,:,c) - b(:,c) ||_2^2  / nbPerCoil <= epsData(c)

              g = @(in) 0.5 * norm( applyWmI( in ), 'fro' )^2;
              gGrad = @(in) applyWmI( applyWmI( in ), 'transp' );
              h = @( in ) indicatorEpsData( applyM( in ) );
              proxth = @(in,t) projectOntoEpsBalls( in );
              if debug == true
                [kFull,objValues,relDiffs] = fista_wLS( kData, g, gGrad, proxth, 'h', h );   %#ok<ASGLU>
              else
                kFull = fista_wLS( kData, g, gGrad, proxth, 'N', nOptIter );
              end
              img = mri_reconRoemer( applyF( kFull, 'inv' ) );

            else  % if numel( epsData ) > 0  &&  max( abs( epsData ) ) ~= 0
              % SPIRiT with hard data consistency
              % minimize (1/2) || W ( Mc^T k + M^T b ) - ( Mc^T k + M^T b ) ||_2^2
              % minimize (1/2) || W Mc^T k + W kData - Mc^T k - kData ||_2^2
              % minimize (1/2) || A k + y ||_2^2 where A = ( W - I ) Mc^T k and y = ( W - I ) kData
              y = applyWmI( kData );
              A_SPIRiT = @(in) applyWmI( applyMc( in, 'transp' ) );
              AH_SPIRiT = @(in) applyMc( applyWmI( in, 'transp' ) );
              g = @(in) 0.5 * norm( A_SPIRiT( in ) + y, 'fro' )^2;
              AHy_SPIRiT = AH_SPIRiT( y );
              gGrad = @(in) AH_SPIRiT( A_SPIRiT( in ) ) + AHy_SPIRiT;

              k0 = zeros( sum( reshape( sampleMask, [], 1 ) == 0 ), 1 );
              % k0 = reshape( k0, [], nCoils );
              % [ys, xs] = ind2sub( sImg, find( sampleMask(:,:,1) == 1 ) );
              % [yq, xq] = ind2sub( sImg, find( sampleMask(:,:,1) == 0 ) );
              % for cIndx = 1 : nCoils
              %   coilKData = kData(:,:,cIndx);
              %   v = coilKData( sampleMask(:,:,cIndx) == 1 );
              %   F = scatteredInterpolant( xs(:), ys(:), v(:), 'linear', 'nearest' );
              %   k0(:,cIndx) = F( xq, yq );
              % end
              % k0( ~isfinite(k0) ) = 0;
              % k0 = k0(:);

              if doChecks == true
                [chk_A,errChk_A] = checkAdjoint( rand( size(k0) ) + 1i * rand( size(k0) ), A, AT );
                if chk_A == true
                  disp( 'Check of A passed' );
                else
                  error([ 'Check of A Adjoint failed with error ', num2str(errChk_A) ]);
                end
              end

              %kFull0 = applyMc( k0(:), 'transp' ) + kData;
              %figure;  showImageCube( 20*log10( abs( kFull0 ) ), 3 );

              if numel( optAlg ) == 0, optAlg = 'lsqr'; end

              if strcmp( optAlg, 'gd' )
                [k,oVals,relDiffs] = gradDescent( k0, gGrad, 'g', g, 'N', nOptIter, 'verbose', verbose, 'printEvery', 50 );   %#ok<ASGLU>
                % Note: SPIRiT requires many iterations; for a sagittal slice of the knee, it required 10,000 iterations.
              elseif strcmp( optAlg, 'lsqr' )
                [ k, lsqrFlag, lsqrRelRes, lsqrIter, lsqrResVec ] = lsqr( ...
                  @applyA_SPIRiT, -y(:), tol, nOptIter, [], [], k0 );   %#ok<ASGLU>
              end
              kFull = reshape( applyMc( k, 'transp' ) + kData, [ sImg nCoils ] );
              img = mri_reconRoemer( applyF( kFull, 'inv' ) );
            end

          else  % if numel( wSize ) > 0

            coilRecons = applyF( kData, 'inv' );
            img = mri_reconRoemer( coilRecons );
          end

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
          img = fista_wLS( img0, g, gGrad, proxth );

        else
          if numel( optAlg ) == 0, optAlg = 'lsqr'; end
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
            img0 = applyP( img0 );
            if nCoils == 1
              [ img, lsqrFlag, lsqrRelRes, lsqrIter, lsqrResVec ] = lsqr( ...
                @applyMFPT, b, tol, nOptIter, [], [], img0(:) );   %#ok<ASGLU>
            else
              [ img, lsqrFlag, lsqrRelRes, lsqrIter, lsqrResVec ] = lsqr( ...
                @apply_vMFPT, b, tol, nOptIter, [], [], img0(:) );   %#ok<ASGLU>
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

function spiritNormWeights = findSpiritNormWeights( kData )

  sampleMask = kData(:,:,1) ~= 0;
  nCoils = size( kData, 3 );

  fftCoords = size2fftCoordinates( size( kData, [1 2] ) );
  ky = fftCoords{1};
  kx = fftCoords{2};
  [ kxs, kys ] = meshgrid( kx, ky );
  kDists = sqrt( kxs.*kxs + kys.*kys );
  sampleMask( kDists == 0 ) = 0;
  kDistsSamples = kDists( sampleMask == 1 );
  [ kDistsSamples, sortedIndxs ] = sort( kDistsSamples );

  spiritNormWeights = cell( 1, 1, nCoils );
  kDistThresh = 0.075;
  parfor coilIndx = 1 : nCoils
    kDataC = kData(:,:,coilIndx);
    F = kDataC( kData(:,:,coilIndx) ~= 0 );
    F = F( sortedIndxs );

    simpleNormWeights = false;
    if simpleNormWeights == true

      coeffs = fitPowerLaw( kDistsSamples, abs( F ), 'ub', [Inf 0], 'linear', false );   %#ok<PFBNS>
      m = coeffs(1);
      p = coeffs(2);
      normWeights = ( kDists.^(-p) ) / m;

    else

      freqBreak = 0.1;

      % Find the parameters that fit the power law for low frequencies    
      coeffs = fitPowerLaw( kDistsSamples( kDistsSamples <= freqBreak ), ...
                            abs( F( kDistsSamples <= freqBreak ) ), 'ub', [Inf 0], 'linear', false );   %#ok<PFBNS>
      mLow = coeffs(1);
      pLow = coeffs(2);
      fitLow = mLow * kDistsSamples.^pLow;
      normWeightsLow = ( kDists.^(-pLow) ) / mLow;

      % Find the parameters tha tfit the power law for high frequencies
      coeffs = fitPowerLaw( kDistsSamples( kDistsSamples >= freqBreak ), ...
                            abs( F( kDistsSamples >= freqBreak ) ), 'ub', [Inf 0], 'linear', false );   %#ok<PFBNS>
      mHigh = coeffs(1);
      pHigh = coeffs(2);
      fitHigh = mHigh * kDistsSamples.^pHigh;
      normWeightsHigh = ( kDists.^(-pHigh) ) / mHigh;

      %figure;  plotnice( kDistsSamples, log( abs( F ) ) );
      %hold on;  plotnice( kDistsSamples, log( fitLow ) );
      %hold on;  plotnice( kDistsSamples, log( fitHigh ) );

      % Blend the parameters together

      %fitBlend = max( fitLow, fitHigh );
      %hold on;  plotnice( kDistsSamples, log( fitBlend ) );
      %figure;  plotnice( kDists(:), normWeights(:) );
      normWeights = min( normWeightsLow, normWeightsHigh );

    end

  
    % fit a line to very small kDists and find the y intercept in order to set the weights when kDists == 0
    smallKDists = kDists( kDists < kDistThresh  &  kDists ~= 0 );
    smallKNormWeights = normWeights( kDists < kDistThresh  &  kDists ~= 0 );
    polyCoeffs = fitPolyToData( 1, smallKDists, smallKNormWeights );
    normWeights( kDists == 0 ) = polyCoeffs(1);

    spiritNormWeights{ coilIndx } = normWeights / sum( normWeights(:) );
  end
  spiritNormWeights = cell2mat( spiritNormWeights );
end


function [w, epsSpiritReg] = findW( acr, wSize, varargin )
  %-- Find the interpolation coefficients w

  p = inputParser;
  p.addOptional( 'spiritNormWeights', [] );
  p.parse( varargin{:} );
  spiritNormWeights = p.Results.spiritNormWeights;

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
        if numel( spiritNormWeights ) > 0
          b( aIndx ) = spiritNormWeights(j,i) * b( aIndx );
          A( aIndx, : ) = spiritNormWeights(j,i) * A( aIndx, : );
        end
        aIndx = aIndx + 1;
      end
    end

    wCoil = A \ b(:);
    epsSpiritReg( coilIndx ) = norm( A * wCoil - b )^2 / numel( b );
    wCoil = [ wCoil( 1 : pt2RemoveIndx - 1 ); 0; wCoil( pt2RemoveIndx : end ); ];
    w{1,1,1,coilIndx} = reshape( wCoil, [ wSize nCoils ] );
  end
  w = cell2mat( w );
end

