
function [ offResMap, phaseOffsetMap ] = mri_mapOffRes( dataCube, TEs, varargin )
  % [ offResMap, phaseOffsetMap ] = mri_mapOffRes( dataCube, TEs [, 'mask', mask ] )
  %
  % Inputs:
  % dataCube - an MxNxK array of data where each k index represents an
  %   image captured at a different echo time.
  % TEs - a 1D array of size K representing the echo times.
  %
  % Optional Inputs:
  % mask - an MxN array.  Any 0 element is data that doesn't get processed.
  %   (default is all 1s)
  %
  % Ouptuts:
  % offResMap - An MxN array representing off resonance frequency (in
  %   radians per unit of echo times)
  % phaseOffsetMap - An MxN array representing the phase at time 0 (in
  %   radians)
  %
  % Written by Nicholas Dwork, Copyright 2018
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addParameter( 'mask', [], @(x) isnumeric(x) || islogical(x) );
  p.parse( varargin{:} );
  mask = p.Results.mask;

  sData = size( dataCube );
  if numel( mask ) == 0, mask=ones( sData(1:2) ); end
  %showImageCube( mask.*abs(dataCube), showScale );  titlenice('t1 ir Data');

  sData = size( dataCube );
  angles = angle( dataCube );

  offResMap = zeros( sData(1:2) );
  phaseOffsetMap = zeros( sData(1:2) );
  for i=1:sData(2)
    for j=1:sData(1)
      if mask(j,i) == 0, continue; end
      unwrapped = squeeze( unwrap( squeeze( angles(j,i,:) ) ) );
      coeffs = fitPolyToData( 1, TEs, unwrapped );
      phaseOffsetMap(j,i) = coeffs(1);
      offResMap(j,i) = coeffs(2);
    end
  end

end

