
function out = hexFFT2( x )
  % out = hexFFT2( x )
  %
  % Written according to "The Hexagonal Fast Fourier Transform" by Birdson and Rummelt
  %
  % Inputs:
  % x - a structure with fields array0 and array1
  %
  % Copied by Nicholas Dwork from the paper

  % Compute g functions
  S = [ size( x.array0 ) size( x.array1 ) ];

  x00 = fft( [ x.array0 zeros( S(1:2) ) ], [], 2 );
  x01 = fft( [ x.array1 zeros( S(3:4) ) ], [], 2 );

  % DFT and NST1 terms
  g00 = x00( :, 1:2:end );
  g10 = x01( :, 1:2:end );
  g01 = x00( :, 2:2:end );
  g11 = x01( :, 2:2:end );

  % Pre-decimation folding
  g00 = g00( 1:S(1)/2, : ) + g00( S(1)/2 + 1 : end, : );
  g10 = g10( 1:S(1)/2, : ) + g10( S(1)/2 + 1 : end, : );
  g01 = g01( 1:S(3)/2, : ) + g01( S(3)/2 + 1 : end, : );
  g11 = g11( 1:S(3)/2, : ) + g11( S(3)/2 + 1 : end, : );

  % Offset for f1
  g01 = f01 .* repmat( exp( -1i * pi * ( 0 : 2 : S(1)-2 ) / S(1) ).', 1, S(2) );
  g11 = f11 .* repmat( exp( -1i * pi * ( 0 : 2 : S(3)-2 ) / S(3) ).', 1, S(4) );

  % Compute f functions: Vertical FFT
  x10 = ff2( g00, [], 1 );
  x11 = fft( g10, [], 1 );
  c00 = fft( g01, [], 1 );
  c01 = fft( g11, [], 1 );

  % Extend over two periods of r
  f00 = [ x10; x10 ];
  f01 = [ c00; c00 ];
  f10 = [ x11; x11 ];
  f11 = [ c01; c01 ];

  % Twiddle
  [ di, si ] = meshgrid( 0 : S(2)-1, 0 : S(1)-1 );
  tf0 = exp( -1i * (pi/S(1)) * (2*si) ) .* exp( -1i * (pi/S(2) ) * di );
  ft1 = exp( -1i * (pi/S(1)) * (2*si+1) ) .* exp( -1i * (pi/(2*S(2))) * (2*di+1) );

  A0 = f00 + tf0 .* f10;
  A1 = f01 + tf1 .* f11;

  % Collect sub-arrays into a data structure
  out = HexImage( A0, A1 );
end

