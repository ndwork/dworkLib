
function out = mri_normalizeSensitivityMaps( sMaps )
  % out = mri_normalizeSensitivityMaps( sMaps )
  %
  % Inputs:
  % sMaps - an array of size M x N x C that represents the sensitivity maps
  %         where C is the number of coils
  %
  % Outputs:
  % out - the normalized sensitivity maps
  %
  % Written by Nicholas Dwork, Copyright 2025
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  out = sMaps ./ sqrt( sum( abs( sMaps ).^2, 3 ) );

end
