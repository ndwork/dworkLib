
function imshownice( img, sdevScale )

  if nargin < 2, sdevScale = 2.5; end;

  meanImg = mean( img(:) );
  sdevImg = std( img(:) );
  
  imshow( img, [ meanImg - sdevScale*sdevImg, meanImg + sdevScale*sdevImg ] );

end
