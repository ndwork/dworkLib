
function mask = mri_makeSampleMask( sMask, nSamples, vdSig, varargin )
  % mask = makeSampleMask( sMask, nSamples, vdSig [, 'maskType', maskType, 'startMask', startMask ] )
  %
  % Inputs:
  % sMask - two element array specifying the size of the mask
  % nSamples - the number of total samples to include in the mask
  % vdSig - standard deviation of distribution as a fraction of sMask
  %
  % Outputs:
  % mask - 2D array of ones and zeros
  %
  % Written by Nicholas Dwork - Copyright 2023
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addParameter( 'maskType', [], @(x) true );
  p.addParameter( 'startMask', zeros( sMask ), @(x) isnumeric(x) || islogical(x) );
  p.parse( varargin{:} );
  maskType = p.Results.maskType;
  startMask = p.Results.startMask;

  if nSamples > prod(sMask), error( 'Too many samples requested.' ); end
  if sum( startMask(:) ) > nSamples, error( 'Too many samples in start mask.' ); end

  mask = startMask;

  for i = 1 : 100

    nSamples2Add = nSamples - sum( mask(:) );
    if nSamples2Add <= 0, break; end

    moreMask = vdSampleMask( sMask, vdSig, nSamples2Add, 'maskType', maskType );

    mask = mask | moreMask;

  end

end



