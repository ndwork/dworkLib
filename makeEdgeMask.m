
function mask = makeEdgeMask( sz, edge )
  % mask = makeEdgeMask( sz, edge )
  %
  % Inputs:
  % sz - two element vector specifying the size of the output image mask
  % edge - scalar specifying number of pixels in edge
  %
  % Outputs:
  % mask - a 2D array where edge pixels are 1 and all other pixels are 0
  %
  % Written by Nicholas Dwork
  %
  % This software is offered under the GNU General Public License 3.0.  It 
  % is offered without any warranty expressed or implied, including the 
  % implied warranties of merchantability or fitness for a particular 
  % purpose.

  if edge < 1, error('Edge must be >= 1'); end;

  mask = zeros( sz );
  mask(1:edge,:) = 1;
  mask(:,1:edge) = 1;
  mask(end-edge+1:end,:) = 1;
  mask(:,end-edge+1:end) = 1;

end

