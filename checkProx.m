
function [ out, err ] = checkProx( x, prox, proxConj, varargin )
  % [ out, err ] = checkProx( x, prox, proxConj [, 'tol', tol, 'nRand', nRand )
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
  sx = size( x );

  v1 = prox( x, 1 );
  v2 = proxConj( x, 1 );
  err = norm( v1(:) + v2(:) - x(:) );
  if err > tol, return; end

  if numel( t ) > 0
    v1 = prox( x, t );
    v2 = t * proxConj( t*x, 1/t );
    err = norm( v1(:) + v2(:) - x(:) );
    if err > tol, return; end
  end

  y = rand( sx );
  v1 = prox( y, 1 );
  v2 = proxConj( y, 1 );
  err = norm( v1(:) + v2(:) - y(:) );
  if err > tol, return; end

  for randIndx = 1 : nRand
    lambda = rand(1,1);
    y = rand( sx );
    v1 = prox( y, lambda );
    v2 = lambda * proxConj( y / lambda, 1 / lambda );
    err = norm( v1(:) + v2(:) - y(:) );
    if err > tol, return; end
  end

  out = true;
end
