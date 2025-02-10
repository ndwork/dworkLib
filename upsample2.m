
function out = upsample2( img, U, varargin )
  % out = upsample2( img, U [, 'S', S, 'sOut', sOut, 'op', op ] )
  %
  % Inputs:
  % img - two dimensional array
  % U - upsample factor.  Either a scalar (assuming same factor in both dimensions)
  %     or an array with two elements (one for each dimension)
  %
  % Optional Inputs:
  % S - amount to shift input, either a scalar or an array with two elements
  % sOut - the size of the output image.  By default, equals size(img) .* U;
  % op - op is either 'notransp' (default) or 'transp', indicating the transpose is desired.
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
  p.addParameter( 'op', [], @(x) true );
  p.addParameter( 'S', zeros( ndims(img), 1 ), @isnumeric );
  p.addParameter( 'sOut', [] );
  p.parse( U, varargin{:} );
  op = p.Results.op;
  S = p.Results.S;
  sOut = p.Results.sOut;

  if isscalar(U), U = U * ones( 2, 1 ); end
  if isscalar(S), S = S * ones( 2, 1 ); end

  sImg = size( img );
  if numel( op ) == 0 || strcmp( op, 'notransp' )
    % upsample
    % Here, img is small and out is large.

    if numel( sOut ) == 0, sOut = sImg(:) .* U(:); end

    yqs = S(1) + (0:sImg(1)-1) * U(1) + 1;
    xqs = S(2) + (0:sImg(2)-1) * U(2) + 1;

    out = zeros( sOut(:)' );
    out( yqs, xqs ) = img;

  elseif strcmp( op, 'transp' )
    % downsample so that operation is transpose
    % Here, img is small and out is large.

    sImg = size( img );
    if numel( sOut ) == 0, sOut = sImg(:) ./ U(:); end

    yqs = S(1) + (0:sOut(1)-1) * U(1) + 1;
    xqs = S(2) + (0:sOut(2)-1) * U(2) + 1;

    out = img( yqs, xqs );

  else
    error( 'Unrecognized operation' );
  end
  
end
