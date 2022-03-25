
function out = downsample2( img, D, varargin )
  % out = downsample2( img, D [, 'S', S ] )
  % Implements the adjoint of upsample2
  %
  % Inputs:
  % img - two dimensional array
  % D - downsample factor.  Either a scalar (assuming same factor in both dimensions)
  %     or an array with two elements (one for each dimension)
  %
  % Optional Inputs:
  % S - amount to shift input, either a scalar or an array with two elements
  %
  % Outputs:
  % out - the downsampled image
  %
  % Written by Nicholas Dwork - Copyright 2020
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 2
    disp( 'Usage:  out = downsample2( img, D [, ''S'', S ] ) ' );
    if nargout > 0, out = []; end
    return;
  end

  p = inputParser;
  p.addRequired( 'D', @isnumeric );
  p.addParameter( 'S', zeros( ndims(img), 1 ), @isnumeric );
  p.parse( D, varargin{:} );
  D = p.Results.D;
  S = p.Results.S;

  if numel(D) == 1, D = D * ones( 2, 1 ); end
  if numel(S) == 1, S = S * ones( 2, 1 ); end

  sImg = size( img );
  yqs = S(1)+1 : D(1)+1 : sImg(1);
  xqs = S(2)+1 : D(2)+1 : sImg(2);

  out = img( yqs, xqs );
end
