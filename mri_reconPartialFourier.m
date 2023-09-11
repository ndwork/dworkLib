
function [ out, phaseOut ] = mri_reconPartialFourier( in, sFSR, varargin )
  % out = mri_reconPartialFourier( in, sFSR [, 'op', op ] )
  %
  % Written according to "Partial k-space Reconstruction" by John Pauly
  %
  % Inputs:
  % in - the input array of size Ny x Nx x nCoils representing the MRI data.
  % sFSR - a scalar specifying the size of the fully sampled region: sFSR x Nx
  %
  % Outputs:
  % 
  %
  % Written by Nicholas Dwork - Copyright 2023
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  Ny = size( in, 1 );
  ys = size2imgCoordinates( Ny );
  centerIndx = find( ys == 0, 1 );

  firstY = centerIndx - ceil( (sFSR-1) / 2 );
  lastY = centerIndx + floor( (sFSR-1) / 2 );
  m = 2 / ( lastY - firstY );
  ramp = m * ys + 1.0;
  ramp( 1 : firstY ) = 0;
  ramp( lastY : end ) = 2;

  inW = bsxfun( @times, in, 2-ramp );
  imgW = fftshift2( ifft2( ifftshift2( inW ) ) );

  inLF = in;  % Low-freq data
  inLF(1:firstY-1,:,:) = 0;
  imgLF = fftshift2( ifft2( ifftshift2( inLF ) ) );
  phaseOut = exp( 1i * angle( imgLF ) );

  img = imgW .* conj( phaseOut );
  out = real( img ) .* phaseOut;
end
