
function [ recon, sMaps ] = mrs_reconRefPeak( kData )
  % [ recon, sMaps ] = mrs_reconRefPeak( kData )
  %
  % Performs the peak reference reconstruction described in "Methodology for
  % improved detection of low concentration metabolites in MRS: optimised 
  % combination of signals from multi-element coil arrays" by Hall et al. (2014)
  %
  % Inputs:
  % kData - array representing the Fourier values of size M x N x nSlices x nCoils x F
  %
  % Written by Nicholas Dwork - Copyright 2021
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  coilSpectrums = fftshift( fftshift( uifft2( kData ), 1 ), 2 );

  freqIndx = ndims( coilSpectrums );
  coilIndx = freqIndx - 1;

  maxImgs = squeeze( max( coilSpectrums, [], freqIndx ) );
  sMaps = bsxfun( @times, maxImgs, 1 ./ sqrt( sum( maxImgs.^2, coilIndx ) ) );

  if ndims( kData ) == 5
    [ M, N, nSlices, nCoils ] = size( sMaps );
    sMapped = bsxfun( @times, reshape( conj( sMaps ), [ M N nSlices nCoils ] ), coilSpectrums );
  else
    [ M, N, nCoils ] = size( sMaps );
    sMapped = bsxfun( @times, reshape( conj( sMaps ), [ M N nCoils ] ), coilSpectrums );
  end

  recon = sum( sMapped, coilIndx );
end


