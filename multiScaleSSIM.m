
function out = multiScaleSSIM( in1, in2, varargin )
  % out = multiScaleSSIM( in1, in2 [, 'k1', k1, 'k2', k2, 'L', L, 'N', N ] )
  %
  % Computes the multi-scale structural similarity metric between inputs 1 and 2 according
  % to "multi-scale structural similarity for image quality assessment" by Wang et al.
  %
  % Inputs:
  % in1 - the first array
  % in2 - the second array of the same size as in1
  %
  % Optional Inputs:
  % k1,k2 - parameters used by the metric
  % L - the dynamic range; if floats and doubles are input, then should be 1
  %     If integers are input, should be 2^# of bits.
  % N - the smallest size considered
  %
  % Written by Nicholas Dwork, Copyright 2021
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 2
    disp( 'Usage:   out = multiScaleSSIM( in1, in2 [, ''k1'', k1, ''k2'', k2, ''L'', L, ''N'', N ] )' );
    if nargout > 0, out = []; end
    return
  end

  p = inputParser;
  p.addParameter( 'k1', 0.01, @isnumeric );
  p.addParameter( 'k2', 0.03, @isnumeric );
  p.addParameter( 'L', [], @ispositive );
  p.addParameter( 'N', 50, @ispositive );
  p.parse( varargin{:} );
  k1 = p.Results.k1;
  k2 = p.Results.k2;
  L = p.Results.L;
  N = p.Results.N;

  dynamicRange = max( [ in1(:); in2(:); ] ) - min( [ in1(:); in2(:); ] );
  if numel( L ) == 0, L = dynamicRange; end

  h = fspecial( 'gaussian', 5, 0.75 );

  sScaled = size( in1 );
  cs = zeros( ceil( log2( min( sScaled ) ) ), 1 );
  ss = zeros( ceil( log2( min( sScaled ) ) ), 1 );

  indxSSIM = 0;
  while min( sScaled ) > N
    indxSSIM = indxSSIM + 1;
    [~,l,c,s] = ssim( in1, in2, 'k1', k1, 'k2', k2, 'L', L );
    cs( indxSSIM ) = c;
    ss( indxSSIM ) = s;

    in1 = imfilter( in1, h, 'circular', 'same' );  in1 = downsample2( in1, 2 );
    in2 = imfilter( in2, h, 'circular', 'same' );  in2 = downsample2( in2, 2 );

    sScaled = floor( sScaled / 2 );
  end

  gamma = 1 / indxSSIM;
  out = l^gamma * prod( c.^gamma .* s.^gamma );
end
