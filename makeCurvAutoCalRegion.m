
function [acr,sLowPass] = makeCurvAutoCalRegion( sImg, varargin )

  p = inputParser;
  p.addOptional( 'beta', 4, @ispositive );
  p.parse( varargin{:} );
  beta = p.Results.beta;

  x = zeros( sImg );
  [~,lowPass] = fdct_wrapping( x, false, 1 );

  sLowPass = size( lowPass );
  kx = kaiser( sLowPass(2), beta );
  ky = kaiser( sLowPass(1), beta );
  [ kxs, kys ] = meshgrid( kx, ky );

  acr = zeros( sImg );
  acr(1:sLowPass(1),1:sLowPass(2)) = kxs .* kys;

  acr = circshift( acr, [ -floor(sLowPass(1) * 0.5) -floor(sLowPass(2) * 0.5) ] );
  acr = fftshift( acr );
end

