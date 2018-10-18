
function R = mri_makeRF( alphas, phis )
  % R = mri_makeRF( alphas, phis )
  % determine the tip-angle rotation matrix applied to (Mx,My,Mz)
  %
  % Inputs:
  % alphas - the tip angles in radians
  %
  % Optional Inputs:
  % phis - the phases of the RF pulse in radians (default is 0)
  %
  % Outputs:
  % the 3x3 rotation matrix
  %
  % NOTE:  All angles in this function follow a left handed convention
  %        I really would have preferred otherwise, but whadareyougonnado?
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 2, phis=0; end;

  calpha = cos(alphas);
  salpha = sin(alphas);
  cphi = cos(phis);  cphiSq = cphi.*cphi;
  sphi = sin(phis);  sphiSq = sphi.*sphi;

  R12 = cphi.*sphi.*(1-calpha);
  R23 = cphi.*salpha;
  R31 = sphi.*salpha;

  if numel( alphas ) == 1

    R = [ ...
     cphiSq + sphiSq*calpha, R12, -R31; ...
     R12, sphiSq+cphiSq*calpha, R23; ...
     R31, -R23, calpha; ...
    ];

  else

    R = zeros( 3, 3, numel(alphas) );
    R(1,1,:) = cphiSq + sphiSq.*calpha;
    R(1,2,:) = R12;
    R(1,3,:) = -R31;
    R(2,1,:) = R12;
    R(2,2,:) = sphiSq + cphiSq.*calpha;
    R(2,3,:) = R23;
    R(3,1,:) = R31;
    R(3,2,:) = -R23;
    R(3,3,:) = calpha;

  end

end
