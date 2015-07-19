
function [gammaBar, gamma] = getGammaH()
  % Returns the gyromagnetic ratio of hydrogen

  gammaBar = 4.257;         % kHz / Gauss
  gamma = 2*pi*gammaBar;    % kRad / Gauss  

end
