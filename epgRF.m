
function Qout = epgRF( Qin, alpha, varargin )
  % Qout = epgRF( Qin, alpha [, phi] )
  % Qin is an 3xN array representing the magnetization state
  %   The rows are F+, F-, and Z
  %   N is the number of Fourier coefficients to store
  % alpha is the tip angle in radians
  % phi is the phase of the rotation axis in radians

  R = epgMakeR( alpha, varargin{:} );

  Qout = R * Qin;

end
