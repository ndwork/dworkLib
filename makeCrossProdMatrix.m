
function out = makeCrossProdMatrix( v )
  % out = makeCrossProdMatrix( v )
  %
  % Inputs:
  % v - a 3 elements array
  %
  % Outputs:
  % out - a 3x3 matrix
  %
  % Written by Nicholas Dwork - Copyright 2016
  
  out = zeros(3);
  out(1,2) = -v(3);
  out(1,3) = v(2);
  out(2,1) = v(3);
  out(2,3) = -v(1);
  out(3,1) = -v(2);
  out(3,2) = v(1);
end
