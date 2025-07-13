
function out = proxL1Complex( in, thresh, weights, b )
  % out = proxL1Complex( in, thresh [, weights, b ] )
  %
  % Returns the proximal operator of f(x) = thresh * L1( x - b )_weights, where
  %   x is a complex vector.
  %
  % Inputs:
  % in - an array of complex values
  % thresh - the thresholding value
  %
  % Optional Inputs:
  % weights - an array of size in specifying threshold scaling factor
  %   for each component (for a weighted L1 norm)
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 2
    disp( 'Usage:  out = proxL1Complex( in, thresh [, weights, b ] )' );
    if nargout > 0, out = []; end
    return;
  end

  if nargin > 2  &&  numel( weights ) > 0
    thresh = thresh .* weights;
  end

  if nargin < 4
    magIn = abs( in );
  else
    magIn = abs( in - b );
  end

  if isreal( in )
    out = sign(in) .* max( magIn - thresh, 0 );

  else
    magOut = max( magIn - thresh, 0 );
    out = zeros( size( in ) );
    out( magIn > thresh ) = in( magIn > thresh ) ./ magIn( magIn > thresh ) .* magOut( magIn > thresh );
  end

  if nargin >= 4
    out = out + b;
  end
end
