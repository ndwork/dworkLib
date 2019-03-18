
function out = smoothImg( in, varargin )
  % out = smoothImg( in [, N, 'gaussian', sigma, 'op', op ] );
  %
  % Inputs:
  % in - either a 2D array (representing an image) or a 3D array (representing a color image)
  %   if op is 'notransp' then in should be the image
  %   if op is 'transp' then adjoint is applied to in
  %
  % Optional Inputs:
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

  nChannels = size( in, 3 );
  out = zeros( size(in) );

  if strcmp( op, 'transp' )
    if mod( N(1), 2 ) == 0, error('N must be odd'); end
    if mod( N(2), 2 ) == 0, error('N must be odd'); end

    for i=1:N(2)
      shiftI = i - floor(N(2)/2) - 1;

      for j=1:N(1)
        shiftJ = j - floor(N(1)/2) - 1;

        for ch=1:nChannels
          shifted = shiftImg( in, [ shiftJ shiftI 0 ] );
          out(:,:,ch) = out + h(j,i) * shifted;
        end
      end
    end

  else
    for ch=1:nChannels
      out(:,:,ch) = imfilter( in(:,:,ch), h, 'same' );
    end

  end
end
