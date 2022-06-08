
function out = evalGaussPDF( x, varargin )
  % out = evalGaussPDF( x [, 'gMean', gMean, 'gVar', gVar, 'gSig', gSig ] )
  % 
  % Evaluation the normal distribution at specific domain values
  %
  % Inputs:
  % x - an array of evaluation points
  %
  % Optional Inputs:
  % gMean - a 1D array or scalar specifying the mean of the distribution
  % gSig - a 1D array or scalar specifying the standard deviation of the distribution
  %   Note that this is slightly slower than using gVar
  % gVar - a 1D array or scalar specifying the variance of the distribution
  %
  % Outputs:
  % out = a 1D array the size of x with values equal to the normal distribution evaluations
  %
  % Written by Nicholas Dwork, Copyright 2019
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addParameter( 'gMean', 0, @isnumeric );
  p.addParameter( 'gSig', [], @isnumeric );
  p.addParameter( 'gVar', 1, @ispositive );
  p.parse( varargin{:} );
  gMean = p.Results.gMean;
  gSig = p.Results.gSig;
  gVar = p.Results.gVar;

  if numel( gSig ) ~= 0, gVar = gSig .* gSig; end

  if gMean ~= 0
    tmp = x - gMean;
    tmp = tmp .* tmp;
  else
    tmp = x .* x;
  end

  out = 1 / ( sqrt( 2 * pi .* gVar ) ) .* exp( -tmp ./ ( 2*gVar ) );
end
