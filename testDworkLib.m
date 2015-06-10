
function testDworkLib

  %% makeDftMatrix
  M = 100;
  out1 = makeDftMatrix( M );
  out2 = makeDftMatrix( M, M );
  diff = out2 - out1;
  err = max( abs( diff(:) ) );
  if err > 1d-12
    error(['makeDftMatrix failed with error ', num2str(err)]);
  else
    disp('makeDftMatrix passed');

end
