
function [ recon, sMaps, lambda ] = mri_reconJointStrucSparseSENSE( kData, varargin )
  % [ recon, sMaps ] = mri_reconJointStrucSparseSENSE( kData [, 'nOuterIter', nOuterIter, 
  %   'noiseCov', noiseCov, 'polyOrder', polyOrder, 'relDiffThresh', relDiffThresh ] );
  %
  % Inputs:
  % kData - a two dimensional array of complex values; uncollected data have values of 0
  %         Its size is [ nKy nKx nCoils ].
  %
  % Optional Inputs:
  % maxOuterIter - a scalar reprenting the maximum number of iterations
  % polyOrder - either a scalar representing the order in both dimensions or
  %             a two element array representing [ yOrder xOrder ]
  % relDiffThresh - dynamic stopping criteria for joint estimation iterations
  %
  % Outputs:
  % recon - a two dimensional complex array that is the reconstructed image
  %
  % Optional Outputs:
  % sMaps - a three dimensional complex array of the sensitivity maps
  %
  % Written by Nicholas Dwork - Copyright 2023
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp([ 'Usage: [ img, sMaps ] = mri_reconJointStrucSparseSENSE( kData [, ''nOuterIter'', nOuterIter, ', ...
      '''noiseCov'', noiseCov, ''reweightEpsilon'', reweightEpsilon, ''nReweightIter'', nReweightIter, ...' ...
      '''polyOrder'', polyOrder,']);
    disp('         ''relDiffThresh'', relDiffThresh, ''verbose'', verbose ] );' )
    if nargout > 0, recon=[]; end
    if nargout > 1, sMaps=[]; end
    return
  end

  p = inputParser;
  p.addParameter( 'lambda', [], @isnumeric );
  p.addParameter( 'noiseCov', [], @isnumeric );
  p.addParameter( 'nReweightIter', 5, @(x) ispositive(x) && mod(x,1)==0 );
  p.addParameter( 'polyOrder', [], @(x) isnonnegative(x) && isinteger(x) );
  p.addParameter( 'relDiffThresh', 0 );
  p.addParameter( 'reweightEpsilon', [], @ispositive );
  p.addParameter( 't', [], @ispositive );
  p.addParameter( 'transformType', 'wavelet', @(x) true );
  p.addParameter( 'waveletType', 'Daubechies-4', @(x) true );
  p.addParameter( 'wavSplit', [], @isnumeric );
  p.addParameter( 'verbose', false, @(x) isnumeric(x) || islogical(x) );
  p.parse( varargin{:} );
  lambda = p.Results.lambda;
  noiseCov = p.Results.noiseCov;
  nReweightIter = p.Results.nReweightIter;
  polyOrder = p.Results.polyOrder;
  relDiffThresh = p.Results.relDiffThresh;
  reweightEpsilon = p.Results.reweightEpsilon;
  t = p.Results.t;
  transformType = p.Results.transformType;
  wavSplit = p.Results.wavSplit;
  waveletType = p.Results.waveletType;
  verbose = p.Results.verbose;

  coilRecons = mri_reconIFFT( kData, 'multiSlice', true );
  roemerRecon = mri_reconRoemer( coilRecons );
  recon = roemerRecon;

  if relDiffThresh > 0
    objValue = [];
    normKData = norm( kData( kData ~= 0 ) );
  end

  for iter = 1 : nReweightIter
    if verbose == true
      disp([ 'Working on Joint Structured Sparsity iteration ', num2str(iter) ]);
    end

    sMaps = mri_makeSensitivityMaps( kData, 'alg', 'sortaSimple' );

    if nargout > 1 && iter == nReweightIter
      [recon,lambda] = mri_reconStructuredSparseSENSE( kData, sMaps, lambda, 'img0', roemerRecon, 'noiseCov', noiseCov, ...
        'optAlg', 'fista_wLS', 'reweightEpsilon', reweightEpsilon, 't', t, ...
        'transformType', transformType, 'waveletType', waveletType, 'wavSplit', wavSplit );
    else
      recon = mri_reconStructuredSparseSENSE( kData, sMaps, lambda, 'img0', roemerRecon, 'noiseCov', noiseCov, ...
        'optAlg', 'fista_wLS', 'reweightEpsilon', reweightEpsilon, 't', t, ...
        'transformType', transformType, 'waveletType', waveletType, 'wavSplit', wavSplit );
    end

    if relDiffThresh > 0
      lastObjValue = objValue;
      sRecons = bsxfun( @times, sMaps, recon );
      fftRecons = fftshift2( fft2( ifftshift2( sRecons ) ) );
      objValue = norm( fftRecons( kData ~= 0 ) - kData( kData ~= 0 ) ) / normKData;

      if iter > 1
        objDiff = ( objValue - lastObjValue ) / lastObjValue;
        if objDiff < relDiffThresh, break; end
      end
    end
  end

end
