
function img = mri_reconModelBased( kData, sMaps )

  coilRecons = mri_reconIFFT( kData );
  img0 = mri_reconRoemer( coilRecons );

  function out = applyS( in, type )
    if nargin < 2 || strcmp( type, 'notransp' )
      out = bsxfun( @times, sMaps, in );
    else
      out = sum( conj(sMaps) .* in, ndims(sMaps) );
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
