
function out = isEven( img, varargin )
  % out = isEven( img, [ threshold ] )
  % Determines whether or not img is even (with circular boundary
  %   conditions)
  % Indexes are defined according to fftshift
  %
  % Inputs:
  % img is a 2D array that is the image
  % threshold is an optional input.  Max magnitude of difference between 
  %   img(x,y) and img(-x,-y) must be less than threshold to be considered
  %   even.

  defaultThresh = 0;
  p = inputParser;
  p.addOptional( 'thresh', defaultThresh, @isnumeric );
  p.parse( varargin{:} );
  thresh = p.Results.thresh;

  mirrorImg = rot90( img, 2 );

  [ny,nx] = size(img);
  nyEven = ~mod(ny,2);
  nxEven = ~mod(nx,2);
  if nyEven || nxEven
    mirrorImg = circshift( mirrorImg, +[nyEven nxEven] );
  end

  absDiff = abs( mirrorImg - img );
  maxDiff = max( absDiff(:) );
  if maxDiff > thresh, out = false; else out = true; end;
end
