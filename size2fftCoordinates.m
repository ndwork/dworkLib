
function k = size2fftCoordinates( N )
  % k = size2fftCoordinates( N )
  %
  % Inputs:
  %   N is an an array specifying the number of elements in each dimension.
  %   For example, N can be a two element array [Ny Nx] specifying the
  %     number of row and columns of the Fourier Transformed image.
  %
  % Outputs:
  %   k(i,:) is a vector with k-space locations for each dimension in
  %     the Fourier domain.
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
  
  firstK = size2fftCoordinates_1D( N(1) );
  k = zeros( numel(firstK), numN );
  k(:,1) = firstK;
  
  for i=2:numN
    k(:,i) = size2fftCoordinates_1D( N(i) );
  end
  
end


function k = size2fftCoordinates_1D( N )
  dk = 1 / N;
  if mod( N, 2 ) == 0
    k = linspace(-0.5,0.5-dk,N);
  else
    k = linspace(-0.5+dk/2,0.5-dk/2,N);
  end
end

