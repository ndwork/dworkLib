
function out = totalVariation( in, varargin )
  % out = totalVariation( in [, lambda ]  )
  %
  % Computes the L2L1 norm (the total variation) of in.
  %
  % Written by Nicholas Dwork - Copyright 2021
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp( 'Usage:  out = totalVariation( in )' );
    if nargout > 0, out=[]; end
    return
  end

  p = inputParser;
  p.addOptional( 'lambda', 1, @positive );
  p.parse( varargin{:} );
  lambda = p.Results.lambda;

  gIn = LpNorms( computeGradient( in ) );
  out = lambda * sum( gIn(:) );
end
