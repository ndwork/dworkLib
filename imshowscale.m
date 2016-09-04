
function imshowscale( img, scale, varargin )
  % imshowscale( img, scale [, method, 'range', range ] )
  % displays figure to the screen where size of image is scaled by scale
  %
  % Inputs:
  % img - 2D array representing the grayscale image or 3D array
  %       representing the color image (third dimension has size 3)
  %
  % Optional Inputs:
  % scale - factor to scale the size of the image for display
  % method - when scaling method, interpolation method to use
  %   default is 'nearest'
  %   any method accepted by imresize is accepted for this parameter
  % range - two element array specifying the display range of intensities
  %   if range is [], sets equal to [min(img(:)) max(img(:))]
  %
  % Written by Nicholas - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  defaultMethod = 'nearest';
  defaultRange = 0;
  p = inputParser;
  p.addOptional( 'method', defaultMethod );
  p.addParameter( 'range', defaultRange );
  p.parse( varargin{:} );
  method = p.Results.method;
  range = p.Results.range;

  if range == 0
    if ismatrix( img )
      imshow( imresize( img, scale, method ) );
    elseif ndims(img) == 3
      imshow( imColorResize( img, scale, method ) );
    else
      error('wrong number of dimensions');
    end
  else
    if ismatrix( img )
      imshow( imresize( img, scale, method ), range );
    elseif ndims(img) == 3
      imshow( imColorResize( img, scale, method ), range );
    else
      error('wrong number of dimensions');
    end
  end
end
