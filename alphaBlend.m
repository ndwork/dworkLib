
function fused = alphaBlend( img1, img2, alpha )
  % out = alphaBlend( img1, img2 [, alpha ] )
  %
  % Perform alpha blending image fusion
  %
  % Inputs:
  % img1 - a 2D or 3D array representing the first image
  % img2 - a 2D or 3D array representing the third image
  %
  % Optional Inputs:
  % alpha - out = alpha * img1 + (1-alpha) * img2
  %
  % Output:
  % fused - the fused image
  %
  % Written by Nicholas Dwork
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if ndims(img1) ~= ndims(img2)
    if ismatrix( img1 )
      img1 = repmat( img1, [1 1 3] );
    elseif ismatrix( img2) 
      img2 = repmat( img2, [1 1 3] );
    end
  end

  if nargin < 3
    alpha = 0.5;
  end

  fused = alpha * img1 + (1-alpha) * img2;
end

