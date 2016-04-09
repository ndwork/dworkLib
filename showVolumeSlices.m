
function showVolumeSlices( vol, varargin )
  % showVolumeSlices( vol, [ range ] );
  %
  % Show the middle slices in each dimension
  %
  % Inputs:
  % vol - a 3D array
  % 
  % Optional Inputs:
  % range - a 2 element array specifying the display range
  %
  % Written by Nicholas Dwork - Copyright 2016

  defaultRange = [];
  p = inputParser;
  p.addOptional( 'range', defaultRange );
  p.parse( varargin{:} );
  displayRange = p.Results.range;

  sVol = size( vol );

  midIndxs = ceil( sVol / 2 );

  slice1 = squeeze( vol( :, :, midIndxs(3) ) );
  slice2 = squeeze( vol( :, midIndxs(2), : ) );
  slice3 = squeeze( vol( midIndxs(1), :, : ) );

  maxY = max( [sVol(1) sVol(2)] );

  padded1 = padData( slice1, [maxY sVol(2)] );
  padded2 = padData( slice2, [maxY sVol(3)] );
  padded3 = padData( slice3, [maxY sVol(3)] );

  out = [ padded1 padded2 padded3 ];
  imshow( out, displayRange );
end
