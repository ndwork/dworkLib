
function [img,relRes] = mri_reconModelBased( kData, sMaps, varargin )
  % [img,relRes] = mri_reconModelBased( kData, sMaps [, 'traj', traj, 'support', support ] )
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
  p.parse( varargin{:} );
  support = p.Results.support;
  traj = p.Results.traj;

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
      out = sum( conj(sMaps) .* in, ndims(sMaps) );
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
      out = applyF( applyS( in ) );
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
      out = applyS( applyF( in, 'transp' ), 'transp' );
      if numel( support ) > 0
        out = out( support ~= 0 );
      end
    end
    out = out(:);
  end

  if numel( traj ) == 0
    coilRecons0 = mri_reconIFFT( kData, 'multiSlice', true );
    kData( dataMask == 1 );
  else
    coilRecons0 = grid_2D( kData, traj, sImg );
  end
  img0 = mri_reconRoemer( coilRecons0 );

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

  function out = g( in )
    Ein = applyE( in );
    if numel( traj ) > 0
      out =  0.5 * real( dotP( Ein - kData(:), Ein - kData(:) ) );
    else
      out =  0.5 * real( dotP( Ein - kData(dataMask == 1), Ein - kData( dataMask == 1 ) ) );
    end
  end

  if numel( support ) > 0
    img0 = img0( support == 1 );
  end

  optAlg = 'lsqr';
  if strcmp( optAlg, 'cgs' )
    [img,optFlag,relres] = cgs( @applyE, kData(:), [], 1000, [], [], [] );   %#ok<ASGLU>
  elseif strcmp( optAlg, 'gradDescent' )
    if numel( traj ) > 0
      ETb = applyE( kData(:), 'transp' );
    else
      ETb = applyE( kData( dataMask == 1 ), 'transp' );
    end
    gGrad = @(x) applyE( applyE( x ), 'transp' ) - ETb;
    [img,objValues,relDiffs] = gradDescent( img0(:), gGrad, 'useLineSearch', true, 'N', 1000, 'g', @g, 'verbose', true );   %#ok<ASGLU>
    %[img,objValues,relDiffs] = gradDescent( img0(:), gGrad, 'useLineSearch', false, 't', 1d-6, 'N', 1000, 'g', @g, 'verbose', true );   %#ok<ASGLU>
  elseif strcmp( optAlg, 'lsqr' )
    [img,optFlag,relRes] = lsqr( @applyE, kData(:), [], 1000, [], [], img0(:) );   %#ok<ASGLU>
  end
  img = reshape( img, [ sImg nSlices ] );
end
