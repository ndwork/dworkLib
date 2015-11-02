
function out = epgMakeR( alpha, phi )
  % alpha is the tip angle (in radians);
  % phi is the phase of the axis at which the RF pulse
  %   is applied (in radians)

  if nargin<2, phi=0; end;

  out = zeros( 3, 3 );

  cosHalfAlpha = cos( alpha/2 );
  cosHalfAlphaSq = cosHalfAlpha * cosHalfAlpha;
  sinAlpha = sin( alpha );
  sinHalfAlpha = sin( alpha/2 );
  sinHalfAlphaSq = sinHalfAlpha * sinHalfAlpha;

  out(1,1) = cosHalfAlphaSq;
  out(2,1) = exp(-1i*2*phi) * sinHalfAlphaSq;
  out(3,1) = -1i/2 * exp(-1i*phi) * sinAlpha;
  out(1,2) = exp(1i*2*phi) * sinHalfAlphaSq;
  out(2,2) = cosHalfAlphaSq;
  out(3,2) = 1i/2 * exp(1i*phi) * sinAlpha;
  out(1,3) = -1i * exp(1i*phi) * sinAlpha;
  out(2,3) = 1i * exp(-1i*phi) * sinAlpha;
  out(3,3) = cos(alpha);
end

