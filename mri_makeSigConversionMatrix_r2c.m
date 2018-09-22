
function M = mri_makeSigConversionMatrix_r2c
  % M = mri_makeSigConversionMatrix_r2c()
  % Create the matrix such that [ M_{xy}; M_{xy}^*; Mz; ] = M [ Mx; My; z; ]
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

  M = [ [ 1  1i 0 ]; ...
        [ 1 -1i 0 ]; ...
        [ 0  0  1 ] ];

end
