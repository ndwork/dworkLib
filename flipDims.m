
function x = flipDims( x, varargin )
  % out = flipDims( in, [, 'dims', dims, 'fftOrigin', true/false ] )
  %
  % Inputs:
  % in - an array of any number of dimensions
  %
  % Optional Inputs:
  % dims - a 1D array specifying the indices to flip.
  %   If not provided, all dimensions are flipped.
  % fftOrigin - if set to true, flips the array about the origin used by fftshift
  %
  % Outputs:
  % out - an array of size equal to in where all dimensions have been flipped
  %
  % Written by Nicholas Dwork - Copyright 2025
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  if nargin < 1
    disp( 'Usage: out = flipDims( in, [, ''dims'', dims, ''fftOrigin'', true/false ] )' );
    if nargout > 0, x = []; end
    return
  end

  p = inputParser;
  p.addParameter( 'dims', [], @ispositive );
  p.addParameter( 'fftOrigin', true );
  p.parse( varargin{:} );
  dims = p.Results.dims;
  fftOrigin = p.Results.fftOrigin;

  if numel( dims ) == 0, dims = 1 : ndims( x ); end

  if fftOrigin == true
    for i = 1 : numel( dims )
      thisDim = dims( i );
      if mod( size( x, thisDim ), 2 ) == 0
        x = circshift( flip( x, thisDim ), 1, thisDim );
      else
        x = flip( x, thisDim );
      end
    end

  else
    for i = 1 : numel( dims )
      x = flip( x, dims(i) );
    end
  end
end
