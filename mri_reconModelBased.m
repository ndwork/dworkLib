
function [img,relRes] = mri_reconModelBased( kData, sMaps, varargin )
  % [img,relRes] = mri_reconModelBased( kData, sMaps [, 'support', support, 'traj', traj ] )
  %
  % img is the argmin of || F S x - b ||_2
  % F is either the FFT or the Non-uniform FFT, based on whether or not traj is supplied
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

  p = inputParser;
  p.addParameter( 'support', [] );
  p.addParameter( 'traj', [], @isnumeric );
  p.addParameter( 'verbose', true );
  p.parse( varargin{:} );
  support = p.Results.support;
  traj = p.Results.traj;
  verbose = p.Results.verbose;

  sImg = size( sMaps, [1 2] );

  sKData = size( kData );
  if numel( traj ) > 0
    nSlices = size( kData, 3 );
  else
    nSlices = size( kData, 4 );
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
      out = applySF( in );
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
      out = applySF( in, 'transp' );
      if numel( support ) > 0
        out = out( support ~= 0 );
      end
    end
    out = out(:);
  end

  img0 = zeros( sImg );

  doCheckAdjoint = false;
  if doCheckAdjoint == true
    innerProd = @(x,y) real( dotP( x, y ) );
    [checkS,errS] = checkAdjoint( img0, @applyS, 'innerProd', innerProd );   %#ok<ASGLU>
    if checkS ~= 1, error( 'Adjoint check of S failed' ); end
    [checkF,errF] = checkAdjoint( repmat(img0,[1,1,8]), @applyF, 'innerProd', innerProd );   %#ok<ASGLU>
    if checkF ~= 1, error( 'Adjoint check of F failed' ); end
    [checkSF,errSF] = checkAdjoint( img0, @applySF, 'innerProd', innerProd );   %#ok<ASGLU>
    if checkSF ~= 1, error( 'Adjoint check of SF failed' ); end
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

  if verbose == true, imgH = figure(); end

  function dispFunc( x, iter )
    if numel( support ) > 0
      tmp = zeros( [ sImg nSlices ] );
      tmp( support ~= 0 ) = x;
      x = tmp;  clear tmp;
    end
    figure( imgH );
    showImageCube( reshape( abs(x), [ sImg nSlices ] ) );
    titlenice([ 'Iteration ', num2str(iter) ]);
  end

  optAlg = 'gradDescent';
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
