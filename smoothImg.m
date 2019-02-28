
function out = smoothImg( in, varargin )
  % out = smoothImg( in [, N, 'gaussian', sigma, 'op', op ] );
  %
  % Inputs:
  % in - if op is 'notransp' then in should be the image
  %      if op is 'transp' then adjoint is applied to in
  % N - the size of each dimension of the smoothing kernel
  % sigma - by default, the smoothing operator is a box car filter.
  %         if sigma is provided, then the smoothing operator is a Gaussian filter
  % op - either 'notransp' (default, meaning filter), or 'transp' (meaning return the
  %      adjoint of the filter)
  %
  % Written by Nicholas - Copyright 2019
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  defaultN = 5;
  defaultSigma = 0;
  p = inputParser;
  p.addOptional( 'N', defaultN );
  p.addParameter( 'gaussian', defaultSigma, @isnumeric );
  p.addParameter( 'op', 'notransp', @(x) true );
  p.parse( varargin{:} );
  N = p.Results.N;
  sigma = p.Results.gaussian;
  op = p.Results.op;

  if numel( N ) == 1, N = [ N N ]; end


  if sigma > 0
    % Gaussian filter with standard deviation sigma
    h = fspecial( 'gaussian', N, sigma );
  else
    % Box car average (or mean) filter
    h = fspecial( 'average', N );
  end

  if strcmp( op, 'transp' )
    if min( mod( N, 2 ) ) == 0, error('smoothImg only works with odd N'); end
    sIn = size( in );
    hN = floor( N / 2 );
    adjIn = zeros( sIn + 2 * hN );

    for i=1:sIn(2)
      for j=1:sIn(1)
        %disp([ j i ]);
        adjIn( j:j+N(1)-1, i:i+N(2)-1 ) = ...
          adjIn( j:j+N(1)-1, i:i+N(2)-1 ) + in(j,i) * h;
      end
    end
    out = adjIn( hN(1)+1 : hN(1)+sIn(1), hN(2)+1 : hN(2)+sIn(2) );

  else
    out = imfilter( in, h );

  end
end
