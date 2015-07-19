
function out = makeDftMatrix( M, N )
  % out = makeDftMatrix( M [, N] )
  %
  % Make the (potentially non-square) DFT matrix for a vector
  % in \mathbb{R}^N

  if nargin<2, N=M; end;

  maxSize = max( M, N );
  fftIn = eye(maxSize);
  fftIn = fftIn(:,1:N);

  out = fft( fftIn );   % fft of each column
end

