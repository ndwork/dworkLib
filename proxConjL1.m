
function out = proxConjL1( in, sigma, t )
  % out = proxConjL1( v, sigma, t )
  %
  % Returns the proximal operator of sigma times the conjugate function of
  % f(x) = t * || x ||_1.
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

  if nargin < 2, sigma = 1; end
  if nargin < 3, t=1; end

  if sigma < 0, error( 'sigma must be non-negative' ); end

  if sigma == 0
    out = in;
    return
  end

  if isreal( in )
    out = min( in, t );
    out = max( out, -t );

  else
    magIn = abs( in );
    scaling = t ./ magIn;
    out = in;
    out( magIn > t ) = in( magIn > t ) .* scaling( magIn > t );

  end

  % out = in - sigma * proxL1Complex( in / sigma, t / sigma );
end
