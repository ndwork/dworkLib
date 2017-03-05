
function Qout = epgGrad( Qin, n )
  % Qin is an 3xN array representing the magnetization state
  %   The rows are F+, F-, and Z
  %   N is the number of Fourier coefficients to store
  % n is 1 / -1 for a positive / negative gradient

  Qout = Qin;
  if abs(n)~=1, error('n must be 1 or -1'); end;

  if n>0
    Qout(1,2:end) = Qin(1,1:end-1);
    Qout(1,1) = conj( Qin(2,2) );
    Qout(2,1:end-1) = Qin(2,2:end);
    Qout(2,end) = 0;

  elseif n<0
    Qout(2,2:end) = Qin(2,1:end-1);
    Qout(2,1) = conj( Qin(1,2) );
    Qout(1,1:end-1) = Qin(1,2:end);
    Qout(1,end) = 0;

  end

end
