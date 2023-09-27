
function out = proxL1Complex( in, thresh, weights, b )
  % out = proxL1Complex( in, thresh [, weights ] )
  %
  % Returns the proximal operator of f(x) = thresh * L1( x - b ), where
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
    disp( 'Usage:  out = proxL1Complex( in, thresh [, weights ] )' );
    return;
  end

  if nargin > 2, thresh = thresh .* weights; end
  if nargin < 4, b = 0; end

  magIn = abs( in - b );
  scalingFactors = thresh ./ magIn;

  out = zeros( size( in ) );
  out( magIn > thresh ) = in( magIn > thresh ) .* ( 1 - scalingFactors( magIn > thresh ) );
  out = out + b;
end
