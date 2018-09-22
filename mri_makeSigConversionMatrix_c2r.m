
function M = mri_makeSigConversionMatrix_c2r
  % M = mri_makeSigConversionMatrix_c2r()
  % Create the matrix such that [ Mx; My; z; ] = M [ M_{xy}; M_{xy}^*; Mz; ]
  %
  % Outputs:
  % M - the 3x3 matrix
  %
  % Written by Nicholas Dwork - Copyright 2018
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  M = [ [  0.5  0.5  0 ]; ...
        [ -0.5i 0.5i 0 ]; ...
        [  0    0    1 ] ];

end
