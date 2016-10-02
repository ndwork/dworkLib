
function imshownice( img, varargin )
  % imshownice( img [, scale, method, 'sdevScale', sdevScale ] )
  % show the image on the following scale
  % meanImg - sdevScale*sdevImg, meanImg + sdevScale*sdevImg
  %
  % Inputs:
  % img - 2D array representing the image
  %
  % Optional Inputs:
  % scale - factor to scale the size of the image for display
  % method - when scaling method, interpolation method to use
  %   default is 'nearest'
  %   any method accepted by imresize is accepted for this parameter
  % sdevScale - default is 2.5
  %
  % Written by Nicholas - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  defaultScale = 1.0;
  defaultMethod = 'nearest';
  defaultSDevScale = 2.5;
  p = inputParser;
  p.addOptional( 'scale', defaultScale );
  p.addOptional( 'method', defaultMethod );
  p.addParameter( 'sdevScale', defaultSDevScale );
  p.parse( varargin{:} );
  scale = p.Results.scale;
  method = p.Results.method;
  sdevScale = p.Results.sdevScale;

  meanImg = mean( img(:) );
  sdevImg = std( img(:) );

  if ismatrix( img )
    % Grayscale image
    tmp = imresize( img, scale, method );
    imshow( tmp, [ meanImg - sdevScale*sdevImg, meanImg + sdevScale*sdevImg ] );
  else
    % Color image
    for i=1:size(img,3)
      tmp(:,:,i) = imresize( img(:,:,i), scale, method );
    end
    inMin = meanImg - sdevScale*sdevImg;
    inMax = meanImg + sdevScale*sdevImg;
    scaled = scaleImg( tmp, [inMin inMax], [0 1] );
    imshow( scaled );
  end
end
