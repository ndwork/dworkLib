
function pcc = phaseCrossCorrelate( img1, img2 )
  % pcc = phaseCrossCorrelate( img1, img2 );
  % Will output the phase cross correlation of image 1 with image 2 as
  % detailed in "An FFT Based Tecnique for Translation, Rotation, and Scale
  % Invariant Image Registration"
  %
  % Inputs:
  % img1 - a 2D array representing an image
  % img2 - a 2D array representing an image
  %
  % Output:
  % pcc - the resulting phase cross correlation
  %
  % Written by Nicholas Dwork - Copyright 2017
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  fftImg1 = fft2( img1 );
  fftImg2 = fft2( img2 );
  numerator = conj(fftImg1) .* fftImg2;
  denominator = abs( numerator );
  fftPC = numerator ./ denominator;
  pcc = ifft2( fftPC );
end
