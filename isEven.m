
function [out,error] = isEven( data, varargin )
  % out = isEven( data, [ threshold ] )
  %
  % Determines whether or not data is even (with circular boundary conditions)
  % Indexes are defined according to fftshift
  %
  % Inputs:
  % data is a 1D or 2D array
  % threshold is an optional input.  The relative error between data(x)
  %   and data(-x) must be less than threshold to be considered even.
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  defaultThresh = 0;
  p = inputParser;
  p.addOptional( 'thresh', defaultThresh, @isnumeric );
  p.parse( varargin{:} );
  thresh = p.Results.thresh;

  numDims = ndims( data );

  if isrow(data) || iscolumn(data)
    [out,error] = isEven1D( data, thresh );

  elseif numDims == 2
    [out,error] = isEven2D( data, thresh );

  end

end


function [out,error] = isEven1D( data, thresh )
  mirrorData = flipud( data(:) );

  nData = numel( data );
  nEven = ~mod( nData, 2 );
  if nEven
    mirrorData = circshift( mirrorData, [nEven 1] );
  end

  error = norm( mirrorData(:) - data(:), 2 ) / norm( data(:), 2 );
  if error > thresh, out = false; else out = true; end   %#ok<SEPEX>
end


function [out,error] = isEven2D( img, thresh )
  [ny,nx] = size(img);
  nyEven = ~mod(ny,2);
  nxEven = ~mod(nx,2);

  mirrorImg = rot90( img, 2 );

  if nyEven || nxEven
    mirrorImg = circshift( mirrorImg, +[nyEven nxEven] );
  end

  error = norm( mirrorImg(:) - img(:), 2 ) / norm( img(:), 2 );
  if error > thresh, out = false; else out = true; end   %#ok<SEPEX>
end

