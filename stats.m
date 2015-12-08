
function stats( data )

  function arrayStats( data )
    minImg = min( data(:) );
    disp(['Min Image: ', num2str(minImg)]);

    maxImg = max( data(:) );
    disp(['Max Image: ', num2str(maxImg)]);

    meanImg = mean( data(:) );
    disp(['Mean Image: ', num2str(meanImg)]);

    stdDevImg = std( data(:) );
    disp(['Std Dev Image: ', num2str(stdDevImg)]);
  end

  imagData = imag(data);
  if max( abs(imagData(:)) ~= 0 )
    disp('Real: ');
    arrayStats( real(data) );
    disp('Imag: ');
    arrayStats( imag(data) );
  else
    arrayStats( data );
  end

end
