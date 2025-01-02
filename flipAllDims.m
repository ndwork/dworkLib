
function x = flipAllDims( x )
  % out = flipAllDims( in )
  %
  % Inputs:
  % in - an array of any number of dimensions
  %
  % Outputs:
  % out - an array of size equal to in where all dimensions have been flipped
  %
  % Written by Nicholas Dwork - Copyright 2025
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  for i = 1 : ndims( x )
    x = flip( x, i );
  end
end
