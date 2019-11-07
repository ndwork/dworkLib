
function out = ifftc( in, varargin )
  % out = ifftc( in [, n, dim ] ) or out = ifftc( in, dim )
  %
  % Compute the centered ifft
  %
  % Inputs:
  % in - an array (of any number of dimensions)
  %
  % Optional Inputs:
  % n - zero pads to this number of elements and then perform an ifft
  %     if a two dimensional ifft, the zero pads to size [n(1) n(2)] and performs ifft
  % dim - if supplied, only does an ifft along this dimension
  % unitary - if true, use a scaling factor to make the transform unitary
  %
  % Written by Nicholas Dwork, Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addOptional( 'n', [], @(x) numel(x) == 0 || ispositive(x) );
  p.addOptional( 'dim', [], @isnumeric );
  p.parse( varargin{:} );
  dim = p.Results.dim;
  n = p.Results.n;

  if numel( dim ) == 0 && numel( n ) ~= 0
    dim = n;  n=[];
  end

  if numel( dim ) > 0
    if numel( n ) > 1, error( 'Too many elements in n with dimension supplied' ); end
    out = ifft( ifftshift( in, dim ), n, dim );

  else

    switch ndims( in )
      case 1
        out = ifft( ifftshift( in ), n );

      case 2
        if numel( n ) > 0
          if numel( n ) ~= 2, error( 'sz has the wrong number of elements' ); end
          out = ifft2( ifftshift( in ), n(1), n(2) );
        else
          out = ifft2( ifftshift( in ) );
        end

      otherwise
        if numel( n ) ~= 0 && numel( n ) ~= ndims( in )
          error( 'n has the wrong number of elements' );
        end
        out = ifftn( ifftshift( in ), n );
    end

  end

end
