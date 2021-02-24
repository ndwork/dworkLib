
function out = shearImg( in, theta, dim )
  % out = shearImg( in, theta [, dim] )
  %
  % shears the image along dimension dim
  %
  % Inputs:
  % in - 2D array
  % theta - rotation angle in radians
  %
  % Optional Inputs:
  % dim - index of dimension to shear
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 3, dim = 1; end

  sIn = size( in );

  midY = ceil((sIn(1)+1)/2);
  midX = ceil((sIn(2)+1)/2);
  ys = (1:sIn(1)) - midY;
  xs = (1:sIn(2)) - midX;
  ys = ys' * ones(1,sIn(2));
  xs = ones(sIn(1),1) * xs;

  if dim == 1
    newYs = xs * tan(-theta) + ys;
    out = interp2( xs, ys, in, xs, newYs, 'linear', 0 );
  else
    newXs = -ys * tan(-theta) + xs;
    out = interp2( xs, ys, in, newXs, ys, 'linear', 0 );
  end

end
