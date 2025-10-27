
function out = undoRadialDistortion( img, ks, varargin )
  p = inputParser;
  p.addParameter( 'c', [] );
  p.addParameter( 'space', [] );
  p.parse( varargin{:} );
  c = p.Results.c;
  space = p.Results.space;
  out = applyRadialDistortion( img, ks, 'c', c, 'dir', -1, 'space', space );
end



