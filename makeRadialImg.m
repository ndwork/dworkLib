
function out = makeRadialImg( sImg, center )
  % makes an image where each point represents the distance from the center
  %
  % out = makeRadialImg( sImg [, center] )
  %
  % Inputs:
  % sImg - 2 element array specifying the (y,x) size of the image
  % center - 2 element array specifying the (y,x) coordinate of the center
  % 
  % Outputs:
  % out - output image
  %
  % Written by Nicholas Dwork - (c) 2016

  y = 1:sImg(1);
  x = 1:sImg(2);

  if nargin < 2
    yCenter = ceil( (sImg(1)+1)/2 );
    xCenter = ceil( (sImg(2)+1)/2 );
  else
    yCenter = center(1);
    xCenter = center(2);
  end
  x = x - xCenter;
  y = y - yCenter;

  [x,y] = meshgrid(x,y);

  out = sqrt( x.*x + y.*y );
end
