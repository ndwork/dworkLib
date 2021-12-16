
function [ recon, sMaps ] = mrs_reconRefPeak( kData, varargin )
  % [ recon, sMaps ] = mrs_reconRefPeak( kData [, 'kTraj', kTraj, 'sImg', sImg ] )
  %
  % Performs the peak reference reconstruction described in "Methodology for
  % improved detection of low concentration metabolites in MRS: optimised 
  % combination of signals from multi-element coil arrays" by Hall et al. (2014)
  %
  % Inputs:
  % kData - array representing the Fourier values of size M x N x nCoils x N1 x N2 x ... x N
  %         or nTraj x nCoils x N1 x N2 x ... x N
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
  % implied warranties of merchantability or fitness for a particular purpose.

  if nargin < 1
    disp( 'Usage:  [ recon, sMaps ] = mrs_reconRefPeak( kData [, ''kTraj'', kTraj ] )' );
    if nargout > 0, recon = []; end
    if nargout > 1, sMaps = []; end
    return
  end

  p = inputParser;
  p.addParameter( 'alpha', [], @isnumeric );
  p.addParameter( 'kTraj', [], @isnumeric );
  p.addParameter( 'nC', [], @ispositive );
  p.addParameter( 'sImg', [], @isnumeric );
  p.addParameter( 'W', [], @ispositive );
  p.parse( varargin{:} );
  alpha = p.Results.alpha;
  kTraj = p.Results.kTraj;
  nC = p.Results.nC;
  sImg = p.Results.sImg;
  W = p.Results.W;

  sKData = size( kData );
  if numel( kTraj ) > 0
    coilSpectrums = mri_reconGridding( kData, kTraj, sImg, 'alpha', alpha, 'W', W, 'nC', nC );
    Ns = sKData(3:end);
  else
    coilSpectrums = fftshift2( uifft2( ifftshift2( kData ) ) );
    Ns = sKData(4:end);
  end
  if numel( Ns ) == 0, Ns = 1; end

  sCSs = size( coilSpectrums );
  coilSpectrums = reshape( coilSpectrums, [ sCSs(1:3) prod( sCSs(4:end) ) ] );

  coilIndx = 3;
  freqIndx = ndims( coilSpectrums );

  maxImgs = zeros( sCSs(1:coilIndx) );
  for u = 1 : size( coilSpectrums, 1 )
    for v = 1 : size( coilSpectrums, 2 )
      for w = 1 : size( coilSpectrums, 3 )
        maxImgs(u,v,w) = max( coilSpectrums( u, v, w, : ), [], freqIndx );
      end
    end
  end

  sMaps = bsxfun( @times, maxImgs, 1 ./ sqrt( sum( abs( maxImgs ).^2, coilIndx ) ) );

  [ M, N, nCoils ] = size( sMaps );
  sMapped = bsxfun( @times, reshape( conj( sMaps ), [ M N nCoils ] ), coilSpectrums );

  recon = sum( sMapped, coilIndx );
  recon = permute( recon, [ 1:coilIndx-1 coilIndx+1:ndims(recon) coilIndx ] );
  sRecon = size( recon );
  recon = reshape( recon, [ sRecon(1:2) Ns ] );
end

