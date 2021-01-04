
function out = ssim( in1, in2, varargin )
  % out = ssim( in1, in2 [, 'k1', k1, 'k2', k2, 'L', L ] )
  %
  % Computes the structural similarity metric between inputs 1 and 2
  %
  % Inputs:
  % in1 - the first array
  % in2 - the second array of the same size as in1
  %
  % Optional Inputs:
  % k1,k2 - parameters used by the metric
  % L - the dynamic range; if floats and doubles are input, then should be 1
  %     If integers are input, should be 2^# of bits.
  %
  % Written by Nicholas Dwork, Copyright 2021
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addParameter( 'k1', 0.01, @isnumeric );
  p.addParameter( 'k2', 0.03, @isnumeric );
  p.addParameter( 'L', 1.0, @ispositive );
  p.parse( varargin{:} );
  k1 = p.Results.k1;
  k2 = p.Results.k2;
  L = p.Results.L;

  mean1 = mean( in1(:) );
  mean2 = mean( in2(:) );
  cov12 = cov( in1(:), in2(:) );
  var1 = cov12(1,1);
  var2 = cov12(2,2);
  cov12 = cov12(1,2);
  c1 = ( k1 * L ).^2;
  c2 = ( k2 * L ).^2;

  num = ( 2 * mean1 * mean2 + c1 ) * ( 2 * cov12 + c2 );
  den = ( mean1*mean1 + mean2*mean2 + c1 ) * ( var1 + var2 + c2 );

  out = num ./ den;
end

