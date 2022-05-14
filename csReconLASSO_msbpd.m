
function [recon,oValues,lambda] = csReconLASSO_msbpd( samples, varargin )
  % recon = csReconLASSO_msbpd( samples, lambda [, 'debug', debug, ...
  %   'nIter', nIter, 'printEvery', printEvery, 'wavSplit', wavSplit, ...
  %   'verbose', verbose, 'waveletType', waveletType ] )
  %
  % This routine minimizes 0.5 * || Ax - b ||_2^2 + lambda || W x ||_1
  %   where A is sampleMask * Fourier Transform * real part, and
  %   W is the wavelet transform.  It assumes a fully sampled center region
  %   corresponding the the two-level sampling scheme defined in
  %   "Breaking the coherence barrier: A new theory for compressed sensing"
  %   by Adcock, Ben, et al.
  %   Moreover, it uses the fully sampled scheme in a way to increase
  %   sparsity in the optimization problem
  %
  % Inputs:
  % samples - a 2D array that is zero wherever a sample wasn't acquired
  % lambda - regularization parameter
  %
  % Optional Inputs:
  % debug - if true, reduces the default number of iterations to 30 and forces verbose
  %         statements during optimization
  % nIter - the number of iterations that FISTA will perform (default is 100)
  % polish - if set to true, adds a polishing step (default is false)
  % printEvery - FISTA prints a verbose statement every printEvery iterations
  % verbose - if true, prints informative statements
  % waveletType - either 'Daubechies' for Daubechies-4 (default) or 'Haar'
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  for vIndx = 1 : numel( varargin )-1
    if strcmp( varargin{vIndx}, 'wavSplit' )
      wavSplit = varargin{ vIndx + 1 };
    end
    if strcmp( varargin{vIndx}, 'transformType' )
      transformType = varargin{ vIndx + 1 };
    end
  end

  sImg = size( samples );
  if strcmp( transformType, 'curvelet' )
    acr = makeCurvAutoCalRegion( sImg );

  elseif strcmp( transformType, 'wavelet' )
    acr = makeWavAutoCalRegion( sImg, wavSplit );

  elseif strcmp( transformType, 'wavCurv' )
    acrCurv = makeCurvAutoCalRegion( sImg );
    acrWav = makeWavAutoCalRegion( sImg, wavSplit );
    acr = max( acrCurv, acrWav );

  else
    error( 'Unrecognized transform type' );
  end

  acrSamplesL = acr .* samples;
  beta = samples - acrSamplesL;
  beta( samples == 0 ) = 0;

  reconL = fftshift( ifft2( ifftshift( acrSamplesL ) ) );

  [reconH,oValues,lambda] = csReconLASSO( beta, varargin{:} );

  recon = reconH + reconL;
end

