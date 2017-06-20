
function out = makeRadialImg( sImg, varargin )
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

  [x,y] = makeCoordinateImgs( sImg, varargin{:} );

  out = sqrt( x.*x + y.*y );
end
