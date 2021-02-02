
function fused = wavFusion( cImg, mImg )
  % fused = wavFusion( cImg, mImg )
  %
  % Inputs:
  % cImg - color or monochrome image
  % mImg - monochrome image
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  split = zeros( 4 );
  split(1,1) = 1;

  sImg = size( cImg );
  nChannels = size( cImg, 3 );
  wtColor = zeros( sImg );
  for i = 1 : nChannels
    wtColor(:,:,i) = wtDaubechies2( cImg(:,:,i), split );
  end

  wtMono = wtDaubechies2( mImg, split );

  fused = zeros( sImg );
  activityMono = abs(wtMono);
  activityColor = abs(wtColor);
  for i = 1 : nChannels
    wtFused = wtColor(:,:,i);
    wtFused( activityMono > activityColor(:,:,i) ) = ...
      wtMono( activityMono > activityColor(:,:,i) );
    fused(:,:,i) = iwtDaubechies2( wtFused, split );
  end

  fused = min( max( fused, 0 ), 1 );
end
