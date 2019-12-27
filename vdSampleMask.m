
function mask = vdSampleMask( sMask, sigmas, nSamples, varargin )
  % mask = vdSampleMask( sMask, sigmas, nSamples [, 'maskType', maskType ] )
  %
  % Creates a variable density sampling mask according to a seperable distribution
  %
  % Inputs:
  % sMask - the size of the mask (a 1D array corresponding to the number
  %   of elements in each dimension)
  %
  % Optional Inputs:
  % maskType - string specifying mask type; can be 'Laplacian' or 'Gaussian'
  %
  % Output:
  % samples - a matrix specifying the sample for each dimension
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp([ 'Usage: mask = vdSampleMask( sMask, sigmas, nSamples ', ...
      '[, ''maskType'', maskType ] )' ]);
    return
  end

  p = inputParser;
  p.addParameter( 'maskType', 'Laplacian', @(x) true );
  p.parse( varargin{:} );
  maskType = p.Results.maskType;

  if numel( sigmas ) == 1, sigmas = sigmas * ones( numel(sMask), 1 ); end

  maskCoordinates = size2imgCoordinates( sMask );

  if iscell( maskCoordinates )
    nDims = numel( maskCoordinates );
  else
    nDims = 1;
  end

  dimSamples = cell( 1, nDims );
  for dimIndx = 1 : nDims
    theseCoordinates = maskCoordinates{ dimIndx };

    if strcmp( maskType, 'Laplacian' )
      samplePDF = evalLaplacePDF( theseCoordinates, 'LSig', sigmas(dimIndx) );
    elseif strcmp( maskType, 'Gaussian' )
      samplePDF = evalGaussPDF( theseCoordinates, 'gSig', sigmas(dimIndx) );
    else
      error( 'Unrecognized mask type' );
    end

    theseSamples = samplesFromPMF( samplePDF, theseCoordinates, nSamples );
    dimSamples{dimIndx} = theseSamples - min( theseCoordinates ) + 1;
  end

  samples1D = sub2ind( sMask, dimSamples{:} );

  mask = zeros( sMask );
  mask( samples1D ) = 1;
end

