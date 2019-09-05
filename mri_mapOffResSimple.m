
function [ offResMap, phaseOffsetMap ] = mri_mapOffResSimple( dataCube, TEs, varargin )
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
  %   radians per unit of echo times) according to a right handed rotation
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

  angleOffset = angle( dataCube(:,:,2) .* conj( dataCube(:,:,1) ) );
  offResMap = angleOffset / ( TEs(2) - TEs(1) );
  offResMap = offResMap .* mask;
  phaseOffsetMap = zeros( size(offResMap) );

end


