
function out = smooth( in, w )

  h = fspecial( 'average', w );
  out = imfilter( in, h );

end
