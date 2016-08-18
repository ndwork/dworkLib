% This software is offered under the GNU General Public License 3.0.  It 
% is offered without any warranty expressed or implied, including the 
% implied warranties of merchantability or fitness for a particular 
% purpose.

function out = makeDftMatrix( M, varargin )
  % Make the (potentially non-square) DFT matrix of size M x N
  %
  % out = makeDftMatrix( M [, N, 'direction', direction] )
  %
  % Inputs:
  % M - number of rows
  % N - number of columns; if not specified then N=M
  %
  % Optional Inputs:
  % direction - if direction == 1, then returns the fft matrix
  %   otherwise, it returns the ifft matrix
  % 
  % Outputs:
  % out - DFT matrix
  %
  % Written by Nicholas Dwork - Copyright 2016
  

  defaultDir = '1';
  p=inputParser;
  p.addOptional( 'N', M );
  p.addParamValue( 'direction', defaultDir );
  p.parse( varargin{:} );
  N = p.Results.N;
  direction = p.Results.direction;

  maxSize = max( M, N );
  
  fftIn = eye(maxSize);
  fftIn = fftIn(:,1:N);

%for i=1:N
% fftIn(:,i) = ifftshift( fftIn(:,i) );
%end
  
  if direction > 0
    out = fft( fftIn );   % fft of each column
  else
    out = ifft( fftIn );   % fft of each column
  end
end

