
function out = mri_findErnstAngle( T1, TR )
  % out = mri_findErnstAngle( T1, TR )
  % The Ernst angle is the angle for maximum signal of a GRE spoiled sequence
  % See Ernst RR, Anderson WA.  Application of Fourier transform spectroscopy to magnetic resonance.
  %   Rev Sci Instrum 1966; 37:93-102.
  %
  % Inputs:
  % T1 and TR must be the same units
  %
  % Outputs:
  % out - the Ernst angle in radians
  %
  % Written by Nicholas Dwork - Copyright 2017
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 2
    disp( 'Usage:  out = mri_findErnstAngle( T1, TR )' );
    return
  end

  out = acos( exp( -TR / T1 ) );  
end


