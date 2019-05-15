
function out = lm( f, fPrime, x0, varargin )
  % out = lm( f, fPrime, x0 [, 'L0', L0, 'nIter', nIter ] )
  %
  % This is an implementation of the Levenberg-Marquardt algorithm that
  % utilizes an analytic expression of the gradient of the objective function.
  %
  % Inputs:
  % f - handle for function to be minimized
  % fPrime - handle for gradient (or derivative) of function to be minimized
  %
  % Optional Inputs:
  % x0 - initial guess
  % L0 - initial value of lambda
  % nMax - Maximum number of iterations
  %
  % Written by Nicholas Dwork

  p = inputParser;
  p.addParameter( 'L0', 1, @isnumeric );
  p.addParameter( 'nIter', 100, @isnumeric );
  p.parse( varargin{:} );
  L0 = p.Results.L0;
  nIter = p.Results.nIter;

  nX = numel( x0 );
  if nX == 1
    out = levenbergMarquardt_var1( f, fPrime, x0, L0, nIter );
  else
    out = levenbergMarquardt_multiVar( f, fPrime, x0, L0, nIter );
  end

end


function out = levenbergMarquardt_var1( f, fPrime, x0, L0, nIter )
  % 1 element variable
  x = x0;
  L = L0;
  for i=1:nIter
    Df = fPrime(x);
    fX = f( x );
    tmp = x - Df / ( lambda + Df * Df ) * f( x );
    fTmp = f( tmp );
    if norm( fTmp, 2 ) < norm( fX, 2 )
      L = 0.8 * L;
      x = tmp;
    else
      L = 2 * L;
    end
  end
  out = x;
end


function out = levenbergMarquardt_multiVar( f, fPrime, x0, L0, nIter )
  % multi-element variable
  x = x0;
  nX = numel( x );
  L = L0;
  for i=1:nIter
    Df = fPrime(x);
    fX = f( x );
    tmp = x - ( Df' * Df + L * eye( nX ) ) \ ( Df' * fX );
    fTmp = f( tmp );
    if norm( fTmp, 2 ) < norm( fX, 2 )
      L = 0.8 * L;
      x = tmp;
    else
      L = 2 * L;
    end
  end
  out = x;
end

