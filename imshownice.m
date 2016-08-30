
function imshownice( img, varargin )
  % imshownice( img [, scale, 'sdevScale', sdevScale ] )
  % show the image on the following scale
  % meanImg - sdevScale*sdevImg, meanImg + sdevScale*sdevImg
  %
  % Written by Nicholas - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  defaultScale = 1.0;
  defaultSDevScale = 2.5;
  p = inputParser;
  p.addOptional( 'scale', defaultScale );
  p.addParameter( 'sdevScale', defaultSDevScale );
  p.parse( varargin{:} );
  sdevScale = p.Results.sdevScale;
  scale = p.Results.scale;

  meanImg = mean( img(:) );
  sdevImg = std( img(:) );

  tmp = imresize( img, scale );
  imshow( tmp, [ meanImg - sdevScale*sdevImg, meanImg + sdevScale*sdevImg ] );
end
