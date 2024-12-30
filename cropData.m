
function out = cropData( data, N )
  % out = cropData( data, N )
  % Crops out the center region of the data.
  % (0,0) is defined according to fftshift
  %
  % Inputs:
  % data - array to be cropped
  % N - specified the size of the cropped image
  %   If N is a scalar, then a cube is extracted
  %   If N is an array of size equal to the number of dimensions of data, 
  %     then the final size is [N(1) N(2) .. N(D)]
  %
  % Written by Nicholas Dwork - Copyright 2017
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It 
  % is offered without any warranty expressed or implied, including the 
  % implied warranties of merchantability or fitness for a particular 
  % purpose.

  if nargin < 2
    disp('Usage: out = cropData( data, N )');
    if nargout > 0, out = []; end
    return
  end

  nD = ndims( data );
  if isscalar(N), N = N * ones(nD,1); end

  sData = size(data);
  subIndxs = cell(nD,1);
  for i=1:nD
    if sData(i) == 1
      subIndxs{i} = 1;
      continue;
    end

    if N(i) > sData(i)
      errorMsg = ['Cropping size ', num2str(N(i)), ' for dimension ', ...
        num2str(i), ' is too large'];
      error(errorMsg);
    end

    halfS = sData(i)/2;
    if mod(sData(i),2)==0
      cy = halfS + 1;
      if mod(N(i),2)==0
        halfN = N(i)/2;
        minIndx = cy - halfN;
        maxIndx = cy + halfN - 1;
      else
        halfN = floor(N(i)/2);
        minIndx = cy - halfN;
        maxIndx = cy + halfN;
      end
    else
      cy = ceil(halfS);
      if mod(N(i),2)==0
        halfN = N(i)/2;
        minIndx = cy - halfN;
        maxIndx = cy + halfN - 1;
      else
        halfN = floor(N(i)/2);
        minIndx = cy - halfN;
        maxIndx = cy + halfN;
      end
    end

    subIndxs{i} = minIndx:maxIndx;
  end

  out = data( subIndxs{:} );
end
