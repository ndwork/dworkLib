
function [ out, err ] = checkProxConj( x, prox, proxConj, varargin )
  % [ out, err ] = checkProxConj( x, prox, proxConj [, 'tol', tol, 'nRand', nRand )
  %
  % Uses the Moreau decomposition to see if the proximal operator and the proximal
  % operation of the conjugate function are correctly related.
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
  p.addParameter( 't', [], @ispositive );
  p.addParameter( 'tol', 1d-6, @ispositive );
  p.parse( varargin{:} );
  nRand = p.Results.nRand;
  t = p.Results.t;
  tol = p.Results.tol;

  out = false;
  err = 0;
  sx = size( x );

  if max( min( x(:) ) ) == 0
    x = ones( sx );
  end

  v1 = prox( x, 1 );
  v2 = proxConj( x, 1 );
  thisErr = norm( v1(:) + v2(:) - x(:) );
  if thisErr > tol, return; end
  err = max( err, thisErr );

  if numel( t ) > 0
    v1 = prox( x, t );
    v2 = t * proxConj( t*x, 1/t );
    thisErr = norm( v1(:) + v2(:) - x(:) );
    if thisErr > tol, return; end
    err = max( err, thisErr );
  end

  y = rand( sx );
  v1 = prox( y, 1 );
  v2 = proxConj( y, 1 );
  thisErr = norm( v1(:) + v2(:) - y(:) );
  if thisErr > tol, return; end
  err = max( err, thisErr );

  for randIndx = 1 : nRand
    lambda = rand(1,1);
    y = rand( sx );
    v1 = prox( y, lambda );
    v2 = lambda * proxConj( y / lambda, 1 / lambda );
    thisErr = norm( v1(:) + v2(:) - y(:) );
    if thisErr > tol, return; end
    err = max( err, thisErr );
  end

  out = true;
end
