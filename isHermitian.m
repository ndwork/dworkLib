
function out = isHermitian( img, varargin )
  % out = isHermitian( img, [ threshold ] )
  % Determines whether or not the real part of img is even and the
  %   imaginary part is odd (with circular boundary conditions)
  % Indexes are defined according to fftshift
  %
  % Inputs:
  % img is a 2D array that is the image
  % threshold is an optional input.  Max magnitude of difference must be 
  %   less than threshold to be considered even.

  imgRealIsEven = isEven( real(img), varargin{:} );
  imgImagIsOdd = isOdd( imag(img), varargin{:} );

  out =  imgRealIsEven && imgImagIsOdd;
end
