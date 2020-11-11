

function [ out, err ] = checkProx( x, prox, f, varargin )
  % [ out, err ] = checkProx( x, prox, f [, 'nRand', nRand )
  %
  % u = prox_f(x)   iff   x - u in subderivative of h at u
  %                 iff   f(z) >= f(u) + (x-u)^T (z-u)
  %
  % Note: if this function fails, then prox is faulty.  If it passes, it does
  %   not guarantee that prox is correct.
  %
  % Inputs:
  % x - argument for the proximal operator
  % prox - function handle to proximal operator that accepts input and scalar
  % f - function handle to the function for which prox is the proximal operator
  %
  % Optional Inputs:
  % nRand - the number of randome vectors to evaluate
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
  p.addParameter( 't', 1, @(x) ispositive(x) || x==0 );
  p.parse( varargin{:} );
  nRand = p.Results.nRand;
  t = p.Results.t;

  out = false;
  err = 0;

  u = prox( x, t );

  for randIndx = 1 : nRand
    z = rand( size( u ) );

    fu = f( u );
    fz = f( z );
    if fz < fu + dotP( z - u, x - u )
      err = fu + dotP( z - u, x - u ) - fz;
      return;
    end
  end

  out = true;
end
