
function out = makeDftMatrix( M, N )

  if nargin < 2
    out = fft(eye(M));
    return;
  end

  out = zeros( M, N );
  j = 0:M-1;
  for i=0:N-1
    out(:,i+1) = exp( -1i * 2*pi * i * j ./ M );
  end

end

