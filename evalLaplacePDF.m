
function out = evalLaplacePDF( x, varargin )
  % out = evalLaplacePDF( x [, 'LMean', LMean, 'LVar', LVar, 'LSig', LSig ] )
  % 
  % Evaluation the Laplace distribution at specific domain values
  %
  % Inputs:
  % x - an array of evaluation points
  %
  % Optional Inputs:
  % LMean - a 1D array or scalar specifying the mean of the distribution
  % LSig - a 1D array or scalar specifyign the standard deviation of the distribution
  %   Note that this is slightly slower than using gVar
  % LVar - a 1D array or scalar specifying the variance of the distribution
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
  p.addParameter( 'LMean', 0, @isnumeric );
  p.addParameter( 'LSig', [], @isnumeric );
  p.addParameter( 'LVar', 1, @ispositive );
  p.parse( varargin{:} );
  LMean = p.Results.LMean;
  LSig = p.Results.LSig;
  LVar = p.Results.LVar;

  if numel( LSig ) ~= 0, LVar = LSig .* LSig; end

  b = sqrt( 0.5 * LVar );

  out = 1 / ( 2.*b ) .* exp( - abs( x - LMean ) / b );

end
