
function [ recon, sMaps ] = mrs_reconRefPeak( kData, varargin )
  % [ recon, sMaps ] = mrs_reconRefPeak( kData [, 'kTraj', kTraj, 'sImg', sImg ] )
  %
  % Performs the peak reference reconstruction described in "Methodology for
  % improved detection of low concentration metabolites in MRS: optimised 
  % combination of signals from multi-element coil arrays" by Hall et al. (2014)
  %
  % Inputs:
  % kData - array representing the Fourier values of size M x N x nSlices x nCoils x N1 x N2 x ... x N
  %         or nTraj x nSlices x nCoils x N1 x N2 x ... x N
  %
  % Optional Inputs:
  % kTraj - An Nx2 array specifying (ky,kx) for each trajectory point
  % sImg - A two element array specifying the size of the image to be reconstructed
  %
  % Written by Nicholas Dwork - Copyright 2021
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp( 'Usage:  [ recon, sMaps ] = mrs_reconRefPeak( kData [, ''kTraj'', kTraj ] )' );
    if nargout > 0, recon = []; end
    if nargout > 1, sMaps = []; end
    return
  end

  p = inputParser;
  p.addParameter( 'kTraj', [], @isnumeric );
  p.addParameter( 'sImg', [], @isnumeric );
  p.parse( varargin{:} );
  kTraj = p.Results.kTraj;
  sImg = p.Results.sImg;

  sKData = size( kData );
  if numel( kTraj ) > 0
    coilSpectrums = mri_reconGridding( kData, kTraj, sImg );
    Ns = sKData(4:end);
  else
    coilSpectrums = fftshift( fftshift( uifft2( kData ), 1 ), 2 );
    Ns = sKData(5:end);
  end

  sCSs = size( coilSpectrums );
  coilSpectrums = reshape( coilSpectrums, [ sCSs(1:4) prod( sCSs(5:end) ) ] );

  coilIndx = 4;
  freqIndx = ndims( coilSpectrums );

  maxImgs = max( coilSpectrums, [], freqIndx );
  sMaps = bsxfun( @times, maxImgs, 1 ./ sqrt( sum( maxImgs.^2, coilIndx ) ) );


  [ M, N, nSlices, nCoils ] = size( sMaps );
  sMapped = bsxfun( @times, reshape( conj( sMaps ), [ M N nSlices nCoils ] ), coilSpectrums );

  recon = sum( sMapped, coilIndx );
  sRecon = size( recon );
  recon = reshape( recon, [ sRecon(1:coilIndx-1) Ns ] );
end

