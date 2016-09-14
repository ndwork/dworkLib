
function pyramid = makeImagePyramid( img, nLevels, spacing )
  % pyramid = makeImagePyramid( img, nLevels, spacing )
  %
  % Inputs:
  % img - 2D array to make a pyramid of
  % nLevels - number of levels of the pyramid
  % spacing - relative size of each pyramid level
  %
  % Outputs:
  % pyramid - a cell array with one element for each level of the pyramid
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 2
    nLevels = 4;
  end
  if nargin < 3
    spacing = 2;
  end

  smooth_sigma = spacing/2;
  f = fspecial('gaussian', 2*round(1.5*smooth_sigma) +1, smooth_sigma);
  ratio = 1. / spacing;

  pyramid = cell(nLevels,1);
  tmp = img;
  for m = 1:nLevels
    pyramid{m} = tmp;
    tmp = imfilter(tmp, f, 'corr', 'symmetric', 'same');  % Gauss filter
    tmp = imresize(tmp, ratio, 'bilinear');  % Downsampling
  end
end