
function out = smoothImg( in, varargin )
  % out = smoothImg( in [, N, 'gaussian', sigma, 'op', op ] );
  %
  % Inputs:
  % in - either a 2D array (representing an image) or a 3D array (representing a color image)
  %   if op is 'notransp' then in should be the image
  %   if op is 'transp' then adjoint is applied to in
  %
  % Optional Inputs:
  % N - either a scalar or a 1D array specifying the size of each dimension of the
  %     smoothing kernel
  % sigma - by default, the smoothing operator is a box car filter.
  %         if sigma is provided, then the smoothing operator is a Gaussian filter
  % op - either 'notransp' (default, meaning filter), or 'transp' (meaning return the
  %      adjoint of the filter)
  %
  % Written by Nicholas - Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp( 'Usage:  out = smoothImg( in [, N, ''gaussian'', sigma, ''op'', op ] );' );
    return
  end

  defaultSigma = 0;
  p = inputParser;
  p.addOptional( 'N', [] );
  p.addParameter( 'gaussian', defaultSigma, @isnumeric );
  p.addParameter( 'op', 'notransp', @(x) true );
  p.parse( varargin{:} );
  N = p.Results.N;
  sigma = p.Results.gaussian;
  op = p.Results.op;

  if numel( N ) == 0
    N = ceil( max( 5*sigma, 3 ) );
    N( mod( N, 2 ) == 0 ) = N( mod( N, 2 ) == 0 ) + 1;
  end

  if min( mod( N, 2 ) ) == 0, error( 'All values of N must be odd' ); end
  if numel( N ) == 1, N = [ N N ]; end

  if sigma > 0
    % Gaussian filter with standard deviation sigma
    h = fspecial( 'gaussian', N, sigma );
  else
    % Box car average (or mean) filter
    h = fspecial( 'average', N );
  end

  nChannels = size( in, 3 );
  out = zeros( size(in) );

  if strcmp( op, 'transp' )

    if mod( N(1), 2 ) == 0
      halfN1 = N(1)/2;
    else
      halfN1 = ceil( N(1)/2 );
    end
    
    if mod( N(2), 2 ) == 0
      halfN2 = N(2)/2;
    else
      halfN2 = ceil( N(2)/2 );
    end

    for i=1:N(2)
      shiftI = i - halfN2;

      for j=1:N(1)
        shiftJ = j - halfN1;

        shifted = shiftImg( in, [ shiftJ shiftI 0 ] );
        out = out + h(j,i) * shifted;
      end
    end

  else

    outCells = cell( 1, 1, nChannels );
    parfor ch = 1 : nChannels
      outCells{ch} = imfilter( in(:,:,ch), h, 'same' );
    end
    out = cell2mat( outCells );

  end
end
