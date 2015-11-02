
function imgStats( img )

  minImg = min( img(:) );
  disp(['Min Image: ', num2str(minImg)]);

  maxImg = max( img(:) );
  disp(['Max Image: ', num2str(maxImg)]);

  meanImg = mean( img(:) );
  disp(['Mean Image: ', num2str(meanImg)]);

  stdDevImg = std( img(:) );
  disp(['Std Dev Image: ', num2str(stdDevImg)]);

end
