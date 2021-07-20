
function stats( data )
  % stats( data )
  %
  % display several relevant stats about the data array to the screen
  %
  % Written by Nicholas Dwork - Copyright 2017
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  function arrayStats( data, preText )
    if nargin < 2, preText = ''; end

    subs = cell(1,ndims(data));

    [minData,minIndx] = min( data(:) );
    [subs{:}] = ind2sub(size(data),minIndx);
    for i=1:numel(subs), subs{i} = num2str(subs{i}); end
    nMinLocs = sum( data(:) == minData );
    if nMinLocs > 1
      disp([ preText, 'Min: ', num2str(minData), '.  There are ', num2str(nMinLocs), ...
        ' locations with this value.' ]);
    else
      minLoc = [ '(', strjoin(subs, ', '), ')' ];
      disp([ preText, 'Min: ', num2str(minData), ' at ', minLoc ]);
    end

    [maxData,maxIndx] = max( data(:) );
    [subs{:}] = ind2sub(size(data),maxIndx);
    for i=1:numel(subs), subs{i} = num2str(subs{i}); end
    nMaxLocs = sum( data(:) == maxData );
    if nMaxLocs > 1
      disp([ preText, 'Max: ', num2str(maxData), '.  There are ', num2str(nMaxLocs), ...
        ' locations with this value.' ]);
    else
      maxLoc = [ '(', strjoin(subs, ', '), ')' ];
      disp([ preText, 'Max: ', num2str(maxData), ' at ', maxLoc ]);
    end

    meanData = mean( data(:) );
    disp([ preText, 'Mean: ', num2str(meanData) ]);

    medData = median( data(:) );
    disp([ preText, 'Median: ', num2str(medData) ]);
    
    stdDevData = std( data(:) );
    disp([ preText, 'Std Dev: ', num2str(stdDevData) ]);

    varData = var( data(:) );
    disp([ preText, 'Var: ', num2str(varData) ]);

    l2norm = norm( data(:), 2 );
    disp([ preText, 'L2 Norm: ', num2str(l2norm)]);
    disp([ preText, 'Mean L2 Norm: ', num2str(l2norm/numel(data)) ]);
  end

  disp([ 'Size of Input: ', num2str(size(data)) ]);
  if numel( data ) == 0, return; end

  imagData = imag(data);
  if max( abs(imagData(:)) ~= 0 )
    disp('Mag: ');
    arrayStats( abs(double(data)), '  ' );
    disp('Phase: ');
    arrayStats( angle(double(data)), '  ' );
    disp('Real: ');
    arrayStats( real(double(data)), '  ' );
    disp('Imag: ');
    arrayStats( imag(double(data)), '  ' );
  else
    arrayStats( double(data) );
  end

end
