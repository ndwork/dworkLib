
function out = randRician( rNu, rSig, varargin )
  % Generate a Rician distribution with parameters rNu and rSig according to
  % https://en.wikipedia.org/wiki/Rice_distribution#Related_distributions
  %
  % out = randRician( rNu, rSig [, theta ] )
  %
  % Inputs:
  % rNu - a single element or an array parameterizing the distribution
  % rSig - a single element or an array (of results size)  specifying the standard deviation of the data.
  %
  % Optional Inputs:
  % theta - a single element or an array (of results size) parameterizing the distribution
  %
  % Written by Nicholas Dwork, 2019

  p = inputParser;
  p.addOptional( 'theta', [], @isnumeric );
  p.parse( varargin{:} );
  theta = p.Results.theta;

  if numel( rNu ) > 1 && numel( rSig ) > 1 && ...
    max( size( rNu ) ~= size( rSig ) ) == 1
    error( 'Inappropriate sizes for rNu and rSig' );
  end

  if numel( rNu ) == 1
    sData = size( rSig );
  else
    sData = size( rNu );
  end

  if numel( theta ) == 0
    theta = pi * rand( sData );
  end

  x = randn( sData ) .* rSig + rNu .* cos(theta);
  y = randn( sData ) .* rSig + rNu .* sin(theta);

  out = sqrt( x.*x + y.*y );
end
