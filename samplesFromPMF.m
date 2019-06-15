
function pmfSamples = samplesFromPMF( pmf, values, varargin )
  % pmfSamples = samplesFromPMF( pmf, values [, N ] )
  %
  % Samples from a specified probability mass function
  %
  % Inputs:
  % values - the image of the random variable  (1D array)
  % pmf - the corresponding probability mass function
  %
  % Optional Inputs:
  % N - the number of samples to generate (default is 1)
  %
  % Outputs:
  % pmfSamples = the sample
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addOptional( 'N', 1, @ispositive );
  p.parse( varargin{:} );
  N = p.Results.N;

  if numel( values ) ~= numel( pmf )
    error( 'values and cdf must have the same number of elements' );
  end

  pmf = pmf(:) / sum( pmf(:) );

  cmf = cumsum( pmf );

  uSamples = rand( N, 1 );  % uniform samples

  pmfSamples = interp1( cmf, values, uSamples, 'nearest', 'extrap' );
end
