
function coords = size2imgCoordinates( N )
  % coords = size2imgCoordinates( N )
  %
  % Inputs:
  %   N is an an array specifying the number of elements in each dimension.
  %   For example, N can be a two element array [Ny Nx] specifying the
  %     number of row and columns of the data.
  %
  % Outputs:
  %   If N is a scalar, then coords is a vector with the coordinates
  %     for each bin of the data.
  %   If N is an array, then coords is a cell array where coords{i} is a
  %     vector with locations for each bin of the i^th dimension.

  if N(1)==1 && numel(N) > 1
    N = N(2:end);
  end
  numN = numel(N);

  coords = cell( numN, 1 );
  for i=1:numN
    coords{i} = size2imgCoordinates_1D( N(i) );
  end

  if numel( coords ) == 1, coords=coords{1}; end
end

function coords = size2imgCoordinates_1D( N )  
  coords = (0:N-1)' - floor(0.5*N);
end
