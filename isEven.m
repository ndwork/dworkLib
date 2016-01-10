
function [out,error] = isEven( img, varargin )
  % out = isEven( data, [ threshold ] )
  % Determines whether or not data is even (with circular boundary
  %   conditions)
  % Indexes are defined according to fftshift
  %
  % Inputs:
  % data is a 1D or 2D array
  % threshold is an optional input.  Relative difference between data(x)
  %   and data(-x) must be less than threshold to be considered even.

  defaultThresh = 0;
  p = inputParser;
  p.addOptional( 'thresh', defaultThresh, @isnumeric );
  p.parse( varargin{:} );
  thresh = p.Results.thresh;

  numDims = ndims( img );
  
  switch numDims
    case 1
      [out,error] = isEven1D( img, thresh );
    case 2
      [out,error] = isEven2D( img, thresh );
  end

end


function [out,error] = isEven1D( data, thresh )
  mirrorData = flipud( data(:) );

  nData = numel( data );

  nEven = -mod( nData, 2 );
  if nEven
    mirrorData = circshift( mirrorData, nEven );
  end

  error = norm( mirrorData(:) - data(:), 2 ) / norm( data(:), 2 );
  if error > thresh, out = false; else out = true; end;
end


function [out,error] = isEven2D( img, thresh )
  [ny,nx] = size(img);
  nyEven = ~mod(ny,2);
  nxEven = ~mod(nx,2);

  mirrorImg = rot90( img, 2 );

  if nyEven || nxEven
    mirrorImg = circshift( mirrorImg, [nyEven nxEven] );
  end

  error = norm( mirrorImg(:) - img(:), 2 ) / norm( img(:), 2 );
  if error > thresh, out = false; else out = true; end;
end

