
function out = fftc( in, varargin )
  % out = fftc( in, n, dim] ) or out = fftc( in, dim] )
  %
  % Compute the centered fft
  %
  % Inputs:
  % in - an array (of any number of dimensions)
  %
  % Optional Inputs:
  % n - zero pads to this number of elements and then perform an fft
  %     if a two dimensional fft, the zero pads to size [n(1) n(2)] and performs fft
  % dim - if supplied, only does an fft along this dimension
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
    out = fftshift( fft( in, n, dim ), dim );

  else

    switch ndims( in )
      case 1
        out = fftshift( fft( in, n ) );

      case 2
        if numel( n ) > 0
          if numel( n ) ~= 2, error( 'sz has the wrong number of elements' ); end
          out = fftshift( fft2( in, n(1), n(2) ) );
        else
          out = fftshift( fft2( in ) );
        end

      otherwise
        if numel( n ) ~= 0 && numel( n ) ~= ndims( in )
          error( 'n has the wrong number of elements' );
        end
        out = fftshift( fftn( in, n ) );
    end

  end

end
