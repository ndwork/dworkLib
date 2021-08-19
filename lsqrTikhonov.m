
function out = lsqrTikhonov( A, b, gamma, varargin )
  % out = lsqrTikhonov( A, b, gamma )
  %
  % Uses lsqr to solve the following optimization problem:
  % minimize || A x - b ||_2^2 + gamma || x ||_2^2
  %
  % Inputs:
  % A - matrix
  % b - vector
  % gamma - scalar regularization paramter
  % All inputs accepted by lsqr are accepted as optional inputs
  %
  % Written by Nicholas Dwork - Copyright 2017
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if gamma<0, error('lsqrTikhonov: gamma must be non-negative'); end

  % M = ( A, I )
  nb = numel(b);
  function out = applyM( in, type )
    if strcmp( type, 'transp' )
      in1 = in(1:nb);
      in2 = in(nb+1:end);
      if isa( A, 'function_handle' )
        ATin1 = A( in1, type );
      else
        ATin1 = A' * in1;
      end
      out = ATin1 + gamma*gamma*in2;
    else
      if isa( A, 'function_handle' )
        Ain = A( in, type );
      else
        Ain = A * in;
      end
      out = [ Ain(:); gamma*gamma*in; ];
    end
  end

  % We will actually minimize || M x - B ||_2 where M=(A,I) and B=(b,0)
  if isa( A, 'function_handle' )
    ATb = A( b, 'transp' );
  else
    ATb = A' * b;
  end
  nX = numel( ATb );
  B = [ b; zeros(nX,1); ];

  out = lsqr( @applyM, B, varargin{:} );
end

