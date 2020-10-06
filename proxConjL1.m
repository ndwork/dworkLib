
function out = proxConjL1( in, t )
  % out = proxConjL1( v, t )
  %
  % Returns the proximal operator of the conjugate function of f(x) = t * || x ||_1.
  % The conjugate function of a norm is the indicator function of the unit ball of
  % the dual norm.  The dual of the L1 norm is the L_Infinity norm.
  % Thus, this operation is a clipping operation that reduces the magnitude of each
  % component to t.
  %
  % Inputs:
  % in - an array of complex values
  % t - the scaling value
  %
  % Outputs:
  % The result of the proximal operator of the conjugate function of the scaled L1 norm
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 2, t=1; end
  if t < 0, error( 't must be positive' ); end

  if t == 0
    out = in;
    return
  end

  tInv = 1 / t;

  if isreal( in )
    out = min( in, tInv );
    out = max( out, -tInv );

  else
    magIn = abs( in );
    scaling = tInv ./ magIn;
    out = in;
    out( magIn > tInv ) = in( magIn > tInv ) .* scaling( magIn > tInv );

  end
end
