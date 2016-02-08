
function coords = size2imgCoordinates( N )
  % coords = size2imgCoordinates( N )
  % Determine the image coordinates where (0,0) is at the "center"
  %   (center is defined in the same way as fftshift)
  %
  % Inputs:
  %   N is an array where the number of elements equal the number of
  %     dimensions of the image.  The value of each element is the size
  %     of the data in that dimension
  %
  % Outputs:
  %   if N is a single element, coords is a 1D array with data coordinates
  %   if N has more than one elements, coords is a cell array where each
  %     element is a 1D array with data coordinates

  numN = numel(N);
  coords = cell(numN,1);
  for i=1:numN
    coords{i} = size2imgCoordinates_1D( N(i) );
  end

  if numel(coords)==1, coords=coords{1}; end;
end


function coords = size2imgCoordinates_1D( N )
  coords = (0:N-1) - floor(0.5*N);
end

