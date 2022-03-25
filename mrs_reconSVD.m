
function out = mrs_reconSVD( coilSpectrums )
  % out = mrs_reconSVD( coilSpectrums )
  %
  % Combines the spectrums from all coils into a single spectral image
  % using the SVD method of "Spectral analysis of multichannel MRS data" by
  % Sandgren et al.
  %
  % Inputs:
  % coilSpectrums - an array of size sImg x nCoils x nF representing the Fourier value
  %   of each spatial location and frequency bin
  %
  % Outputs:
  % out - an array of size sImg x nF
  %
  % Written by Nicholas Dwork, Copyright 2022
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  sCoilSpectrums = size( coilSpectrums );
  ndims = numel( sCoilSpectrums );
  sImg = sCoilSpectrums( 1 : ndims - 2 );
  nCoils = size( coilSpectrums, ndims-1 );
  nF = size( coilSpectrums, ndims );
  nVoxels = prod( sImg );

  coilSpectrums = reshape( coilSpectrums, [ nVoxels nCoils nF ] );

  out = zeros( [ nVoxels nF ] );
  for vox = 1 : nVoxels
    theseSpectrums = transpose( squeeze( coilSpectrums( vox, :, : ) ) );
    [u,s,~] = svd( theseSpectrums, 'econ', 'vector' );
    out( vox, : ) = s( 1 ) * u( : , 1 );
  end

  out = reshape( out, [ sImg nF ] );
end
