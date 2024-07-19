
function mask = vdpdSampleMask( sMask, nSamples, varargin )
  % mask = vdpdSampleMask( sMask, nSamples [, 'Delta', Delta, 'maskType', maskType, 'startMask', startMask ] )
  %
  % Creates a variable density Poisson Disc sampling mask according to a seperable distribution
  %
  % Inputs:
  % sMask - the size of the mask (a 1D array corresponding to the number
  %   of elements in each dimension)
  % sigmas - a 1D array or scalar specifying the standard deviations of the distribution
  %   in each dimension
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
    disp( 'Usage: mask = vdpdSampleMask( sMask, nSamples )' );
    return
  end

  p = inputParser;
  p.addParameter( 'Delta', [], @ispositive );
  p.addParameter( 'maxIter', 100, @ispositive );
  p.addParameter( 'startMask', zeros( sMask ) );
  p.addParameter( 'verbose', true );
  p.parse( varargin{:} );
  Delta = p.Results.Delta;
  maxIter = p.Results.maxIter;
  startMask = p.Results.startMask;
  verbose = p.Results.verbose;

  if numel( Delta ) == 0, Delta = 0.3; end

  if sum( startMask(:) ~= 0 ) >= nSamples, error( 'startMask more than nSamples' ); end

  nDiffSamples = @(in) nSamples - sum( sum( makeSampleMask( in, sMask, startMask, Delta ) ) );

  m = binarySearch( nDiffSamples, 20, 500, 'tol', 0.1, 'verbose', verbose );
  mask = makeSampleMask( m, sMask, startMask, Delta );

  if sum( mask(:) ) > nSamples
    newMaskIndxs = find( ( mask - startMask ) ~= 0 );
    newMaskIndxs = newMaskIndxs( randperm( numel(newMaskIndxs) ) );
    nStartMask = sum( startMask(:) );
    n2add = nSamples - nStartMask;
    mask = startMask;
    mask( newMaskIndxs( 1 : n2add ) ) = 1;
  end

  for sampleIter = 1 : maxIter

    nSamples2Add = nSamples - sum( mask(:) );
    if nSamples2Add <= 0, break; end

    iterMask = makeSampleMask( m, sMask, startMask, Delta );

    uMask = mask | iterMask;  % union of masks

    if sum( sum( mask | iterMask ) ) > nSamples
      % Need to keep just the number needed to get to nSamples
      newSampleIndxs = find( ( uMask - mask)  == 1 );
      nNewSamples = numel( newSampleIndxs );
      newSampleIndxs = newSampleIndxs( randperm( nNewSamples ) );
      iterMask(:) = 0;
      iterMask( newSampleIndxs( 1 : nSamples2Add ) ) = 1;
    end

    mask = mask | iterMask;

    if sum( mask(:) ) == nSamples, break; end
  end

end


function out = makeSampleMask( m, sMask, startMask, Delta )

  dks = 1 ./ ( sMask - 1 );

  %r = @(x) max( ( ( LpNorms( x, 2, 2 ) + 0.2 ) / 100 ), min_r );
  r = @(x) ( LpNorms( x, 2, 2 ) + Delta ) / m;
  min_r = Delta / m;

  ks = poissonDisc2( r, 'min_r', min_r )';
  [ ~, iterSamples ] = movePointsToGrid( bsxfun( @minus, ks, 0.5*dks' ), [-0.5, -0.5], 0.5-dks, sMask );

  out = ( iterSamples > 0 ) | startMask;
end

