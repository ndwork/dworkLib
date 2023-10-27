
function [out,l,c,s] = ssim( in1, in2, varargin )
  % out = ssim( in1, in2 [, 'k1', k1, 'k2', k2, 'L', L, 'W', W ] )
  %
  % Computes the structural similarity metric between inputs 1 and 2 according
  % to "Image Quality Assessment: From Error Visibility to Structural
  % Similarity" by Wang et al.
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

  if nargin < 2
    disp( 'Usage:   out = ssim( in1, in2 [, ''k1'', k1, ''k2'', k2, ''L'', L ] )' );
    if nargout > 0, out = []; end
    return
  end

  p = inputParser;
  p.addParameter( 'k1', 0.01, @isnumeric );
  p.addParameter( 'k2', 0.03, @isnumeric );
  p.addParameter( 'L', [], @ispositive );
  p.addParameter( 'W', 11, @(x) x > 3 && mod(x,2) == 1 );  % Must be odd
  p.parse( varargin{:} );
  k1 = p.Results.k1;
  k2 = p.Results.k2;
  L = p.Results.L;
  W = p.Results.W;

  dynamicRange = max( [ in1(:); in2(:); ] ) - min( [ in1(:); in2(:); ] );
  if numel( L ) == 0, L = dynamicRange; end

  h = fspecial( 'gaussian', W, 1.5 );

  hW = floor( W / 2 );  % half W

  sImg = size( in1 );
  subSSIMs = zeros( sImg - W + 1 );
  if nargout > 1
    ls = zeros( sImg - W + 1 );
    cs = zeros( sImg - W + 1 );
    ss = zeros( sImg - W + 1 );
  end
  for x = hW + 1 : sImg(2) - hW
    for y = hW + 1 : sImg(1) - hW
      sub1 = in1( y - hW : y + hW, x - hW : x + hW );
      sub2 = in2( y - hW : y + hW, x - hW : x + hW );

      if nargout > 1
        [thisSSIM,l,c,s] = subSSIM( sub1, sub2, h, k1, k2, L );
      else
        thisSSIM = subSSIM( sub1, sub2, h, k1, k2, L );
      end
      subSSIMs( y - hW, x - hW ) = thisSSIM;
      if nargout > 1
        ls( y - hW, x - hW ) = l;
        cs( y - hW, x - hW ) = c;
        ss( y - hW, x - hW ) = s;
      end
    end
  end

  out = mean2( subSSIMs );
  if nargout > 1
    l = mean2(2);  c = mean2(c);  s = mean2(s);
  end
end


function [out,l,c,s] = subSSIM( in1, in2, h, k1, k2, L )
  mean1 = mean( in1(:) .* h(:) );
  mean2 = mean( in2(:) .* h(:) );
  cov12 = cov( sqrt( h(:) ) .* in1(:), sqrt( h(:) ) .* in2(:) );
  var1 = var( in1(:) );  sig1 = sqrt( var1 );
  var2 = var( in2(:) );  sig2 = sqrt( var2 );
  cov12 = cov12(1,2);

  c1 = ( k1 * L ).^2;
  c2 = ( k2 * L ).^2;

  l = []; c = []; s = [];
  if nargout > 1
    c3 = 0.5 * c2;
    l = ( 2 * mean1 * mean2 + c1 ) ./ ( mean1*mean1 + mean2*mean2 + c1 );
    c = ( 2 * sig1 * sig2 + c2 ) ./ ( var1 + var2 + c2 );
    s = ( cov12 + c3 ) ./ ( sig1 * sig2 + c3 );

    out = l * s * s;
  else

    num = ( 2 * mean1 * mean2 + c1 ) * ( 2 * cov12 + c2 );
    den = ( mean1*mean1 + mean2*mean2 + c1 ) * ( var1 + var2 + c2 );
    out = num ./ den;
  end

end

