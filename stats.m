
function stats( data )

  function arrayStats( data, preText )
    if nargin < 2, preText = ''; end;

    subs = cell(1,ndims(data));

    [minData,minIndx] = min( data(:) );
    [subs{:}] = ind2sub(size(data),minIndx);
    for i=1:numel(subs), subs{i} = num2str(subs{i}); end;
    minLoc = [ '(', strjoin(subs, ', '), ')' ];
    disp([preText, 'Min: ', num2str(minData), ' at ', minLoc ]);

    [maxData,maxIndx] = max( data(:) );
    [subs{:}] = ind2sub(size(data),maxIndx);
    for i=1:numel(subs), subs{i} = num2str(subs{i}); end;
    maxLoc = [ '(', strjoin(subs, ', '), ')' ];
    disp([preText, 'Max: ', num2str(maxData), ' at ', maxLoc ]);

    meanData = mean( data(:) );
    disp([preText, 'Mean: ', num2str(meanData)]);

    stdDevData = std( data(:) );
    disp([preText, 'Std Dev: ', num2str(stdDevData)]);
    
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
