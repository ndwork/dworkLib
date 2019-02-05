
function out = randRician( rNu, rSig )
  % Generate a Rician distribution with parameters rNu and rSig according to
  % https://en.wikipedia.org/wiki/Rice_distribution#Related_distributions
  %
  % Inputs:
  % rNu - a single element or an array parameterizing the distribution
  % rSig - a single element or an array specifying the standard deviation of the data.
  %   Note that if rNu and rSig are both arrays, they must have the same size.
  %
  % Written by Nicholas Dwork, 2019

  if numel( rNu ) > 1 && numel( rSig ) > 1 && ...
    max( size( rNu ) ~= size( rSig ) ) == 1
    error( 'Inappropriate sizes for rNu and rSig' );
  end

  if numel( rNu ) == 1
    sData = size( rSig );
  else
    sData = size( rNu );
  end

  theta = pi * rand( sData );

  x = randn( sData ) .* rSig + rNu .* cos(theta);
  y = randn( sData ) .* rSig + rNu .* sin(theta);

  out = sqrt( x.*x + y.*y );
end
