
function out = calcAutofocusMetric( in, varargin )
  % Calculates the autofocus metric according to "Blind retrospective motion correction of
  % MR images" by Loktyushin et al.
  %
  % Input:
  % in - an array of any number of dimensions
  %
  % Output:
  % out - the value of the autofocus metric
  %
  % Written by Nicholas Dwork, Copyright 2022
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  p = inputParser;
  p.addParameter( 'minDenom', 1d-6, @ispositive );
  p.addParameter( 'minV', 1d-6, @ispositive );
  p.parse( varargin{:} );
  minDenom = p.Results.minDenom;
  minV = p.Results.minV;

  gradIn = computeGradient( in );

  nDimsIn = ndims( in );
  gradIn = reshape( gradIn, [ numel( in ) nDimsIn ] );

  denoms = transpose( diag( gradIn' * gradIn ) ) + minDenom;

  v = sqrt( bsxfun( @rdivide, gradIn .* conj( gradIn ), denoms ) );

  out = sum( diag( -v' * log( v + minV ) ) );
end

