
function out = makeDftMatrix( M, N )

  if nargin<2, N=M; end;

  maxSize = max( M, N );
  fftIn = eye(maxSize);
  fftIn = fftIn(:,1:N);

  out = fft( fftIn );
end

