
function [img,relRes] = mri_reconModelBased( kData, varargin )
  % [img,relRes] = mri_reconModelBased( kData [, sMaps, 'sImg', sImg, 'support', support, 'traj', traj ] )
  %
  % img is the argmin of || F S x - b ||_2 when support is not supplied
  % F is either the FFT or the Non-uniform FFT, based on whether or not traj is supplied
  % When support is supplied, img is the argmin of || F S D^T x - b ||_2
  % where D is the transformation that isolates the values of the image in the support
  %
  % Inputs:
  % For one slice,
  %   kData - Either 1) a 3D array of size M x N x nCoils representing the kSpace data
  %     or 2) an nTraj x nCoils array represent M kSpace data values
  %   sMaps - a 3D array of size M x N x nCoils representing the sensitivity maps
  % For multiple slices
  %   kData - Either 1) (Cartesian) a 3D array of size M x N x nCoils x S representing the kSpace data
  %     or 2) (non-Cartesian) an nTraj x nCoils x S array represent M kSpace data values
  %   sMaps - a 3D array of size M x N x nCoils x S representing the sensitivity maps
  %
  % Optional Inputs:
  % support - a 2D array of 1s and 0s indicating the support of the img
  % traj - by default, assumes Cartesian sampling.  If traj is supplied, then
  %   non-Cartesian sampling is assumed.
  %
  % Written by Nicholas Dwork, Copyright 2023
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp([ 'Usage:  [img,relRes] = mri_reconModelBased( kData [, sMaps, ''sImg'', sImg, ', ...
           '''support'', support, ''traj'', traj ] )' ]);
    if nargout > 0, img = []; end
    if nargout > 1, relRes = []; end
    return
  end

  p = inputParser;
  p.addOptional( 'sMaps', [], @(x) isnumeric(x) || numel(x) == 0 );
  p.addParameter( 'doCheckAdjoint', false );
  p.addParameter( 'optAlg', 'lsqr', @(x) true );
  p.addParameter( 'sImg', [], @(x) numel(x) == 0  ||  numel(x) == 2 );
  p.addParameter( 'showScale', 3 );
  p.addParameter( 'support', [] );
  p.addParameter( 'traj', [], @isnumeric );
  p.addParameter( 'verbose', true );
  p.parse( varargin{:} );
  sMaps = p.Results.sMaps;
  doCheckAdjoint = p.Results.doCheckAdjoint;
  optAlg = p.Results.optAlg;
  sImg = p.Results.sImg;
  showScale = p.Results.showScale;
  support = p.Results.support;
  traj = p.Results.traj;
  verbose = p.Results.verbose;

  if numel( sMaps ) == 0
    if numel( traj ) == 0
      sImg = size( kData, [1 2] );
    elseif numel( sImg ) == 0
      if numel( support ) > 0
        sImg = size( support );
      else
        error( 'Either sImg or support must be supplied for single-coil non-Cartesian reconstruction' );
      end
    end
  else
    sImg = size( sMaps, [1 2] );
  end

  sKData = size( kData );
  if numel( traj ) > 0
    nCoils = size( kData, 2 );
    nSlices = size( kData, 3 );
  else
    nCoils = size( kData, 3 );
    nSlices = size( kData, 4 );
  end

  if nCoils > 1  &&  numel( sMaps ) == 0
    error( 'You must pass in sMaps when nCoils > 1' );
  end

  if ismatrix( support ) && nSlices > 1
    % The support regions for all slices are the same
    support = repmat( support, [ 1 1 nSlices ]);
  end

  function out = applyS( in, type )
    if nargin < 2 || strcmp( type, 'notransp' )
      out = bsxfun( @times, sMaps, in );
    else
      out = sum( conj(sMaps) .* in, 3 );
    end
  end

  function out = applyF( in, type )
    if nargin < 2 || strcmp( type, 'notransp' )
      if numel( traj ) > 0
        out = iGrid_2D( in, traj );
      else
        out = fftshift2( fft2( ifftshift2( in ) ) );
      end
    else
      if numel( traj ) > 0
        out = iGridT_2D( in, traj, sImg );
      else
        out = fftshift2( fft2h( ifftshift2( in ) ) );
      end
    end
  end

  function out = applySF( in, type )
    if nargin < 2 || strcmp( type, 'notransp' )
      out = applyF( applyS( in ) );
    else
      out = applyS( applyF( in, 'transp' ), 'transp' );
    end
  end

  if numel( traj ) == 0, dataMask = ( abs(kData) ~= 0 ); end

  function out = applyE( in, type )
    if nargin < 2 || strcmp( type, 'notransp' )
      if numel( support ) > 0
        tmp = zeros( [ sImg nSlices ] );
        tmp( support ~= 0 ) = in;
        in = tmp;  clear tmp;
      else
        in = reshape( in, [ sImg nSlices ] );
      end
      if numel( sMaps ) > 0
        out = applySF( in );
      else
        out = applyF( in );
      end
      if numel( traj ) == 0
        out = out( dataMask == 1 );
      end
    else
      if numel( traj ) == 0
        tmp = zeros( sKData );
        tmp( dataMask == 1 ) = in;
        in = tmp;  clear tmp;
      end
      in = reshape( in, sKData );
      if numel( sMaps ) > 0
        out = applySF( in, 'transp' );
      else
        out = applyF( in, 'transp' );
      end
      if numel( support ) > 0
        out = out( support ~= 0 );
      end
    end
    out = out(:);
  end

  img0 = zeros( sImg );

  if doCheckAdjoint == true
    innerProd = @(x,y) real( dotP( x, y ) );
    [checkF,errF] = checkAdjoint( repmat(img0,[1,1,nCoils]), @applyF, 'innerProd', innerProd );   %#ok<ASGLU>
    if checkF ~= 1, error([ 'Adjoint check of F failed with error, ', num2str(errF) ]); end
    if numel( sMaps ) > 0
      [checkS,errS] = checkAdjoint( img0, @applyS, 'innerProd', innerProd );   %#ok<ASGLU>
      if checkS ~= 1, error( 'Adjoint check of S failed' ); end
      [checkSF,errSF] = checkAdjoint( img0, @applySF, 'innerProd', innerProd );   %#ok<ASGLU>
      if checkSF ~= 1, error( 'Adjoint check of SF failed' ); end
    end
    if numel( support ) > 0
      [checkE,errE] = checkAdjoint( img0( support == 1 ), @applyE, 'innerProd', innerProd );   %#ok<ASGLU>
    else
      [checkE,errE] = checkAdjoint( img0, @applyE, 'innerProd', innerProd );   %#ok<ASGLU>
    end
    if checkE ~= 1, error( 'Adjoint check of E failed' ); end
  end

  if numel( traj ) > 0
    b = kData(:);
  else
    b = kData( dataMask == 1 );
  end

  function out = g( in )
    EinMb = applyE( in ) - b;
    out =  0.5 * real( dotP( EinMb, EinMb ) );
  end

  if numel( support ) > 0, img0 = img0( support == 1 ); end

  if verbose == true  &&  strcmp( optAlg, 'gradDescent' ), imgH = figure(); end

  function dispFunc( x, iter )
    if numel( support ) > 0
      tmp = zeros( [ sImg nSlices ] );
      tmp( support ~= 0 ) = x;
      x = tmp;  clear tmp;
    end
    figure( imgH );
    showImageCube( reshape( abs(x), [ sImg nSlices ] ), showScale );
    titlenice([ 'Iteration ', num2str(iter) ]);
  end

  if strcmp( optAlg, 'cgs' )
    [img,optFlag,relres] = cgs( @applyE, b, [], 1000, [], [], [] );   %#ok<ASGLU>
  elseif strcmp( optAlg, 'gradDescent' )
    ETb = applyE( b, 'transp' );
    gGrad = @(x) applyE( applyE( x ), 'transp' ) - ETb;
    [img,objValues,relDiffs] = gradDescent( img0(:), gGrad, 'useLineSearch', true, 'N', 100, 'g', @g, ...
      'dispFunc', @dispFunc, 'verbose', verbose );   %#ok<ASGLU>
  elseif strcmp( optAlg, 'lsqr' )
    [img,optFlag,relRes] = lsqr( @applyE, b, [], 250, [], [], img0(:) );   %#ok<ASGLU>
  end

  if numel( support ) > 0
    tmp = zeros( [ sImg nSlices ] );
    tmp( support ~= 0 ) = img;
    img = tmp;  clear tmp;
  else
    img = reshape( img, [ sImg nSlices ] );
  end
end
