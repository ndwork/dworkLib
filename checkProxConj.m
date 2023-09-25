
function [ out, err ] = checkProxConj( x, prox, proxConj, varargin )
  % [ out, err ] = checkProxConj( x, prox, proxConj [, 'tol', tol, 'nRand', nRand )
  %
  % Uses the Moreau decomposition to see if the proximal operator and the proximal
  % operation of the conjugate function are correctly related.
  %
  % Assumes the proximal operator has the form prox( x, t );
  % Assumes the proximal operator of the conjugate function has the form prox( x, sigma, t )
  %
  % Written by Nicholas Dwork, Copyright 2020
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addParameter( 'nRand', 1, @ispositive );
  p.addParameter( 'sigma', 1, @isnonnegative );
  p.addParameter( 't', [], @ispositive );
  p.addParameter( 'tol', 1d-6, @ispositive );
  p.parse( varargin{:} );
  nRand = p.Results.nRand;
  sigma = p.Results.sigma;
  t = p.Results.t;
  tol = p.Results.tol;

  out = false;
  err = 0;
  sx = size( x );

  if max( min( x(:) ) ) == 0
    x = ones( sx );
  end

  v1 = sigma * prox( x/sigma, 1/sigma );
  v2 = proxConj( x, sigma, 1 );
  thisErr = norm( x(:) - v1(:) - v2(:) ) / norm( x );
  err = max( err, thisErr );
  if thisErr > tol, return; end

  if numel( t ) > 0 && t ~= 1
    v1 = sigma * prox( x/sigma, t/sigma );
    v2 = proxConj( x, sigma, t );
    thisErr = norm( x(:) - v1(:) - v2(:) ) / norm( x );
    err = max( err, thisErr );
    if thisErr > tol, return; end
  end

  y = rand( sx );
  v1 = sigma * prox( y/sigma, 1/sigma );
  v2 = proxConj( y, sigma, 1 );
  thisErr = norm( y(:) - v1(:) - v2(:) ) / norm( x );
  err = max( err, thisErr );
  if thisErr > tol, return; end

  for randIndx = 1 : nRand
    randSigma = rand(1,1);
    lambda = rand(1,1);
    y = rand( sx );
    v1 = randSigma * prox( y/randSigma, lambda/randSigma );
    v2 = proxConj( y, randSigma, lambda );
    thisErr = norm( y(:) - v1(:) - v2(:) );
    err = max( err, thisErr );
    if thisErr > tol, return; end
  end

  out = true;
end
