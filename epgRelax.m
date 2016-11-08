function Qout = epgRelax( Qin, t, T1, T2, M0 )
  % Qout = epgRelax( Qin, t, T1, T2, [ M0 ] )
  % Qin is an 3xN array representing the magnetization state
  %   The rows are F+, F-, and Z
  %   N is the number of Fourier coefficients to store
  if nargin < 5, M0=1; end;

  E1 = exp(-t/T1);
  Qout = [ exp(-t/T2) * Qin(1:2,:); ...
           E1 * Qin(3,:); ];
  Qout(3,1) = Qout(3,1) + M0*(1-E1);
end
