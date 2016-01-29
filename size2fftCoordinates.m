
function ks = size2fftCoordinates( N )
  % ks = size2fftCoordinates( N )
  %
  % Inputs:
  %   N is an an array specifying the number of elements in each dimension.
  %   For example, N can be a two element array [Ny Nx] specifying the
  %     number of row and columns of the Fourier Transformed image.
  %
  % Outputs:
  %   If N is a scalar, then ks is a vector with the k-space locations
  %     for each bin in the Fourier domain.
  %   If N is an array, then ks is a cell array where ks{i} is a vector
  %     with k-space locations for each bin of the i^th dimension.
  %
  %   For our 2D example, 
  %     let fftImg = fftshift( fft2( ifftshift( img ) ) );
  %     ky is a vector of size( img, 1 ) with the k-space location of
  %       each column of fftImg
  %     kx is a vector of size( img, 2 ) with the k-space location of
  %       each row of fftImg


  if N(1)==1 && numel(N) > 1
    N = N(2:end);
  end
  numN = numel(N);

  ks = cell( numN, 1 );
  for i=1:numN
    ks{i} = size2fftCoordinates_1D( N(i) );
  end

  if numel( ks ) == 1, ks=ks{1}; end;
end


function k = size2fftCoordinates_1D( N )
  dk = 1 / N;
  if mod( N, 2 ) == 0
    k = transpose( linspace(-0.5,0.5-dk,N) );
  else
    k = transpose( linspace(-0.5+dk/2,0.5-dk/2,N) );
  end
end

