
function out = dotP( a, b )
  % out = dotP( a, b );
  % Computes the dot product between two vectors a and b.
  % a . b = sum( a .* conj(b) )
  %
  % Written by Nicholas Dwork

  out = sum( a(:) .* conj(b(:)) );

end
