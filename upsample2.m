
function out = upsample2( img, U, varargin )
  % out = upsampleImg( img, U [, 'S', S, 'sOut', sOut ] )
  %
  % Inputs:
  % img - two dimensional array
  % U - upsample factor.  Either a scalar (assuming same factor in both dimensions)
  %     or an array with two elements (one for each dimension)
  %
  % Optional Inputs:
  % S - amount to shift input, either a scalar or an array with two elements
  % sOut - the size of the output image.  By default, equals size(img) .* U;
  %
  % Outputs:
  % out - the upsampled image of size sOut
  %
  % Written by Nicholas Dwork - Copyright 2020
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addRequired( 'U', @isnumeric );
  p.addParameter( 'S', zeros( ndims(img), 1 ), @isnumeric );
  p.addParameter( 'sOut', [] );
  p.parse( U, varargin{:} );
  S = p.Results.S;
  sOut = p.Results.sOut;

  sImg = size( img );

  if numel(U) == 1, U = U * ones( 2, 1 ); end
  if numel(S) == 1, S = S * ones( 2, 1 ); end
  if numel( sOut ) == 0, sOut = sImg(:) .* U(:); end

  out = zeros( sOut(:)' );
  yqs = S(1) + (0:sImg(1)-1) * (U(1)+1) + 1;
  xqs = S(2) + (0:sImg(2)-1) * (U(2)+1) + 1;
  out( yqs, xqs ) = img;
end
