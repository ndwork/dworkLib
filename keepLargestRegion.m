
function out = keepLargestRegion( in )
  % out = keepLargestRegion( in )
  %
  % Isolates the largest continuous region
  %
  % Inputs:
  % in - a binary 2D array with values of 1 or 0
  %
  % Outputs:
  % out - a binary 2D array with values of 1 or 0

  [ labels, nLabels ] = bwlabel( in );
  nPerLabel = zeros( nLabels, 1 );
  for labelIndx = 1 : nLabels
    nPerLabel( labelIndx ) = ( sum( find( labels(:) == labelIndx ) ) );
  end

  [~,largestLabel] = max( nPerLabel );
  out = zeros( size( in ) );
  out( labels == largestLabel ) = 1;
end
