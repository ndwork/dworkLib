
function R = mri_makeRF( alpha, phi )

  if nargin < 2, phi=0; end;

  calpha = cos(alpha);
  salpha = sin(alpha);
  cphi = cos(phi);  cphiSq = cphi*cphi;
  sphi = sin(phi);  sphiSq = sphi*sphi;

  R = [ ...
    cphiSq + sphiSq*calpha, cphi*sphi*(1-calpha), -sphi*salpha; ...
    cphi*sphi*(1-calpha), sphiSq+cphiSq*calpha, cphi*salpha; ...
    sphi*salpha, -salpha*cphi, calpha; ...
  ];

end
