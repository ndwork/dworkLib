
function R = mri_makeRF( alpha, phi )
  % R = mri_makeRF( alpha, phi )
  % determine the tip-angle rotation matrix applied to (Mx,My,Mz)
  %
  % Inputs:
  % alpha - the tip angle in radians
  % phi - the phase of the RF pulse in radians
  %
  % Outputs:
  % the 3x3 rotation matrix
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 2, phi=0; end;

  calpha = cos(alpha);
  salpha = sin(alpha);
  cphi = cos(phi);  cphiSq = cphi*cphi;
  sphi = sin(phi);  sphiSq = sphi*sphi;

  R12 = cphi*sphi*(1-calpha);
  R23 = cphi*salpha;
  R31 = sphi*salpha;

  R = [ ...
    cphiSq + sphiSq*calpha, R12, -R31; ...
    R12, sphiSq+cphiSq*calpha, R23; ...
    R31, -R23, calpha; ...
  ];

end
