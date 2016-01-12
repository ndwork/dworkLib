
function [out,error] = isHermitian( data, varargin )
  % [out,error] = isHermitian( data, [ threshold ] )
  % Determines whether or not the real part of data is even and the
  %   imaginary part is odd (with circular boundary conditions)
  % Indexes are defined according to fftshift
  %
  % Inputs:
  % data is a 1D or 2D array
  % threshold is an optional input.  The relative error
  %   must be less than threshold to be considered Hermitian.

  defaultThresh = 0;
  p = inputParser;
  p.addOptional( 'thresh', defaultThresh, @isnumeric );
  p.parse( varargin{:} );
  thresh = p.Results.thresh;

  numDims = ndims( data );

  if isrow(data) || iscolumn(data)
    [out,error] = isHermitian1D( data, thresh );

  elseif numDims == 2
    [out,error] = isHermitian2D( data, thresh );

  end

end


function [out,error] = isHermitian1D( data, thresh )
  mirrorData = flipud( data(:) );

  nData = numel( data );

  nEven = ~mod(nData,2);
  if nEven
    mirrorData = circshift( mirrorData, [nEven 1] );
  end

  diff = data - conj(mirrorData);
  error = norm(diff(:),2) / norm(data(:),2);
  if error > thresh, out = false; else out = true; end;
end


function [out,error] = isHermitian2D( img, thresh )
  mirrorImg = rot90( img, 2 );

  [ny,nx] = size(img);
  nyEven = ~mod(ny,2);
  nxEven = ~mod(nx,2);
  if nyEven || nxEven
    mirrorImg = circshift( mirrorImg, +[nyEven nxEven] );
  end

  diff = img - conj(mirrorImg);
  error = norm(diff(:),2) / norm(img(:),2);
  if error > thresh, out = false; else out = true; end;
end

