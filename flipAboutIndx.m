
function out = flipAboutIndx( in, indx )
  % out = flipAboutIndx( in, indx )
  %
  % Flips the array in about the location of indx
  %
  % Inputs:
  % in - array (may be multi-dimensional)
  % indx - a 1D array specifying the index for each dimension of in
  %
  % Outputs:
  % out - an array of the size of in with its elements flipped
  %
  % Written by Nicholas Dwork - Copyright 2021
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  shifted = circshift( in, -( indx - 1 ) );
  
  flipped = shifted;
  for dim = 1 : ndims( in )
    flipped = flip( flipped, dim );
  end

  out = circshift( flipped, indx ); 
end
