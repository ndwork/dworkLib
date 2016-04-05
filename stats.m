
function stats( data )

  function arrayStats( data, preText )
    if nargin < 2, preText = ''; end;

    minImg = min( data(:) );
    disp([preText, 'Min Image: ', num2str(minImg)]);

    maxImg = max( data(:) );
    disp([preText, 'Max Image: ', num2str(maxImg)]);

    meanImg = mean( data(:) );
    disp([preText, 'Mean Image: ', num2str(meanImg)]);

    stdDevImg = std( data(:) );
    disp([preText, 'Std Dev Image: ', num2str(stdDevImg)]);
    
    l2norm = norm( data(:), 2 );
    disp([preText, 'L2 Norm: ', num2str(l2norm)]);
  end

  imagData = imag(data);
  if max( abs(imagData(:)) ~= 0 )
    disp('Real: ');
    arrayStats( real(double(data)), '  ' );
    disp('Imag: ');
    arrayStats( imag(double(data)), '  ' );
  else
    arrayStats( double(data) );
  end

end
