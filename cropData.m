
function out = cropData( data, N, varargin )
  % out = cropData( data, N [, C ] )
  % Crops out the center region of the data.
  % (0,0) is defined according to fftshift
  %
  % Inputs:
  % data - array to be cropped
  % N - specified the size of the cropped array
  %   If N is a scalar, then a cube is extracted
  %   If N is an array of size equal to the number of dimensions of data, 
  %     then the final size is [N(1) N(2) .. N(D)]
  %
  % Optional Inputs:
  % C - specify the center of the cropped region (same format as N)
  %
  % Output:
  % out - the cropped array
  %
  % Written by Nicholas Dwork - Copyright 2017
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It 
  % is offered without any warranty expressed or implied, including the 
  % implied warranties of merchantability or fitness for a particular 
  % purpose.

  if nargin < 2
    disp('Usage: out = cropData( data, N [, ''C'', C ] )');
    if nargout > 0, out = []; end
    return
  end

  p = inputParser;
  p.addParameter( 'C', [], @isnumeric );
  p.parse( varargin{:} );
  C = p.Results.C;
  clear p;

  nD = ndims( data );
  if isscalar(N), N = N * ones(nD,1); end
  if isscalar(C), C = C * ones(nD,1); end

  sData = size( data );
  if numel( C ) == 0
    halfS = sData / 2;
    C = ceil( halfS + 1 );
  end

  subIndxs = cell(nD,1);
  for i=1:nD
    if sData(i) == 1
      subIndxs{i} = 1;
      continue;
    end

    if N(i) > sData(i)
      error([ 'Cropping size ', num2str(N(i)), ' for dimension ', num2str(i), ' is too large' ]);
    end

    halfN = floor( N(i) / 2 );
    minIndx = C(i) - halfN;
    maxIndx = C(i) + halfN;
    if mod( N(i), 2 ) == 0
      maxIndx = maxIndx - 1;
    end

    if minIndx < 1  ||  maxIndx > sData(i)
      error([ 'Cropping outside of array in dimension ', num2str(i) ]);
    end

    subIndxs{i} = minIndx : maxIndx;
  end

  out = data( subIndxs{:} );
end
