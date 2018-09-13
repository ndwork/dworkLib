
function [gammaBar, gamma] = getGammaH()
  % Returns the gyromagnetic ratio of hydrogen
  % [gammaBar, gamma] = getGammaH()
  %   gammaBar is in kHz / Gauss
  %   gamma is in kRad / Gauss

  gammaBar = 4.2576;        % kHz / Gauss
  gamma = 2*pi*gammaBar;    % kRad / Gauss / s 

end
