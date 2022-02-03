
function out = retinexLCE( img, varargin )
  % out = retinexLCE( img )
  %
  % Performs local contrast enhancement according to "Retinex Processing for Automatic Image
  % Enhancement" by Rahman, Jobson, and Woodell in Journal of Electronic Imaging in 2004.
  %
  % Inputs:
  % img - a 2D array representing an image.  (Note that the paper describes a color version, but
  %       that is not yet implemented)
  %
  % Outputs:
  % out - a 2D array that is the image with local contrast enhancement
  %
  % Written by Nicholas Dwork - Copyright 2022
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    if nargout > 0, out = []; end
    disp( 'Usage:  out = retinexLCE( img )' );
    return;
  end

  p = inputParser;
  p.addParameter( 'sigmas', [ 15, 80, 250 ], @ispositive );
  p.parse( varargin{:} );
  sigmas = p.Results.sigmas;

  logImg = log( img );
  fftImg = fftshift( fft2( ifftshift( img ) ) );

  imgCoords = size2imgCoordinates( size( img ) );
  yCoords = imgCoords{1} * ones( 1, size( img, 2 ) );
  xCoords = ones( size( img, 1 ), 1 ) * imgCoords{2}';
  kerCoords = -( xCoords .* xCoords + yCoords .* yCoords );

  out = zeros( size( img ) );
  nSigmas = numel( sigmas );
  for k = 1 : nSigmas
    sigmaSq = sigmas( k ) * sigmas( k );

    ker = exp( kerCoords / sigmaSq );
    ker = ker / sum( ker(:) );
    fftKer = fftshift( fft2( ifftshift( ker ) ) );

    kerConvI = fftshift( ifft2( ifftshift( fftKer .* fftImg ) ) );

    thisOut = logImg - log( kerConvI );
    out = out + thisOut;
  end

  out = exp( out / nSigmas );
end

