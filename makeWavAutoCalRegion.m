
function [acr,sACR] = makeWavAutoCalRegion( sImg, wavSplit, varargin )

  p = inputParser;
  p.addOptional( 'beta', 4, @ispositive );
  p.parse( varargin{:} );
  beta = p.Results.beta;

  wavMask = makeLowFreqWavMask( sImg, wavSplit );

  lastX = find( wavMask(1,:), 1, 'last' );
  lastY = find( wavMask(:,1), 1, 'last' );

  kx = kaiser( lastX, beta );
  ky = kaiser( lastY, beta );
  [ kxs, kys ] = meshgrid( kx, ky );

  acr = zeros( size( wavMask ) );
  acr(1:lastY,1:lastX) = kxs .* kys;
  %acr(1:lastY,1:lastX) = 1;

  acr = circshift( acr, [ -floor(lastY * 0.5) -floor(lastX * 0.5) ] );
  acr = fftshift( acr );

  if nargout > 1
    sACR = sImg ./ ( 2.^( log( size(wavSplit) * 2 ) ./ log(2) ) );
  end
end


