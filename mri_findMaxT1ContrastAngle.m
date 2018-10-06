
function out = mri_findMaxT1ContrastAngle( T1, TR )
  % out = mri_findMaxT1ContrastAngle( T1, TR )
  %
  % This function calculates the angle that yields maximum contrast
  % It is from the paper entitled, "Optimization of Flip Angle for 
  % T1 Dependent Contrast in MRI"
  %
  % Inputs:
  % T1 and TR must be the same units
  %
  % Outputs:
  % out - the angle of maximum contrast (in radians)
  %
  % Written by Nicholas Dwork - Copyright 2018
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 2
    disp( 'Usage:  out = mri_findMaxT1ContrastAngle( T1, TR )' );
    return
  end

  E1 = exp( -TR / T1 );
  out = acos( ( 2*E1 - 1 ) / ( 2 - E1 ) );
end


