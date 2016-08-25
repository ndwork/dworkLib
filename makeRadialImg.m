
function out = makeRadialImg( sImg )
  % makes an image where each point represents the distance from the center
  %
  % out = makeRadialImg( sImg )
  %
  % Inputs:
  % sImg - 2 element array specifying the size of the image
  % 
  % Outputs:
  % out - output image
  %
  % Written by Nicholas Dwork - (c) 2016

  y = 1:sImg(1);
  x = 1:sImg(2);

  x = x - ceil( (sImg(2)+1)/2 );
  y = y - ceil( (sImg(1)+1)/2 );

  [x,y] = meshgrid(x,y);

  out = sqrt( x.*x + y.*y );
end
