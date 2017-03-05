
function out = epgMakeR( alpha, phi )
  % out = epgMakeR( alpha, phi )
  % The function creates a matrix that multiplies (Mxy,Mxy^*,Mz) to create a left
  % handed precession about the axis designated by angle phi.
  %
  % Inputs:
  % alpha is the tip angle (in radians);
  % phi is the phase of the axis at which the RF pulse is applied (in radians)
  %
  % Outputs:
  % out is a 3x3 rotation matrix applied to EPG coefficients or
  %   (Mxy,conj(Mxy),Mz)
  %
  % Written by Nicholas Dwork - Copyright 2017
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.
  

  if nargin<2, phi=0; end;

  out = zeros( 3, 3 );

  cosHalfAlpha = cos( alpha/2 );
  cosHalfAlphaSq = cosHalfAlpha * cosHalfAlpha;
  sinAlpha = sin( alpha );
  sinHalfAlpha = sin( alpha/2 );
  sinHalfAlphaSq = sinHalfAlpha * sinHalfAlpha;

  out(1,1) = cosHalfAlphaSq;
  out(2,1) = exp(-1i*2*phi) * sinHalfAlphaSq;
  out(3,1) = 0.5i * exp(-1i*phi) * sinAlpha;
  out(1,2) = conj( out(2,1) );
  out(2,2) = cosHalfAlphaSq;
  out(3,2) = -0.5i * exp(1i*phi) * sinAlpha;
  out(1,3) = 1i * exp(1i*phi) * sinAlpha;
  out(2,3) = -1i * exp(-1i*phi) * sinAlpha;
  out(3,3) = cos(alpha);
end

