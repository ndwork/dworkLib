
function out = showVolumeSlices( vol, varargin )
  % showVolumeSlices( vol, [ range, 'rots', rots, 'noshow', noshow ] );
  %
  % Show the middle slices in each dimension
  %
  % Inputs:
  % vol - a 3D array
  % rots - a 3 element array specifying the rotation to apply to each frame
  % noshow - if true, no image is displayed
  % 
  % Optional Inputs:
  % range - a 2 element array specifying the display range
  %
  % Written by Nicholas Dwork - Copyright 2016

  defaultRange = [];
  defaultRots = [ pi/2 pi/2 pi/2 ];
  p = inputParser;
  p.addOptional( 'range', defaultRange );
  p.addParamValue( 'noshow', false );
  p.addParamValue( 'rots', defaultRots );
  p.parse( varargin{:} );
  displayRange = p.Results.range;
  rots = p.Results.rots;
  noshow = p.Results.noshow;

  sVol = size( vol );

  midIndxs = ceil( sVol / 2 );

  slice1 = squeeze( vol( :, :, midIndxs(3) ) );
  slice2 = squeeze( vol( :, midIndxs(2), : ) );
  slice3 = squeeze( vol( midIndxs(1), :, : ) );

  slice1 = imrotate( slice1, rots(1), 'loose' );
  slice2 = imrotate( slice2, rots(2), 'loose' );
  slice3 = imrotate( slice3, rots(3), 'loose' );

  s1 = size( slice1 );
  s2 = size( slice2 );
  s3 = size( slice3 );
  maxY = max([ s1(1) s2(1) s3(1) ]);

  padded1 = padData( slice1, [maxY s1(2)] );
  padded2 = padData( slice2, [maxY s2(2)] );
  padded3 = padData( slice3, [maxY s3(2)] );

  out = [ padded1 padded2 padded3 ];
  if ~noshow
    imshow( out, displayRange );
  end
end
