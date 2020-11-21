
function cube2Gif( cube, varargin )
  % cube2Gif( cube, [ outFile, delayTime, 'dim', dim, 'loopCount', loopCount ] )
  %
  % This function converts a 3D array of data into a gif.  The third dimension
  % is converted into the temporal dimension.
  %
  % Inputs: 
  % cube - a 3D array.  The dynamic range of cube must be [0,1].  The data will be
  %   scaled by 255 before writing to file as 8-bit values.
  % outFile (optional) - a string contianing the name of the gif file
  %   By default, outFiles is assumed ot be 'makeAGifOut.gif'
  % delayTime (optional) - the time between frames in the animation
  %   By default, delayTime is assumed to be 0.3 seconds
  % loopCount (optional) - the number of times to loop the animation
  %   By default, loopCount is infinity (meaning an endless loop)
  %
  % Written by Nicholas Dwork - Copyright 2017
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp( 'Usage:  cube2Gif( cube, [ outFile, delayTime, ' );
    disp( '  ''dim'', dim, ''loopCount'', loopCount ] )' );
    return
  end
  
  defaultOutFile = 'cube2GifOut.gif';
  defaultDelayTime = 0.3;
  defaultLoopCount = Inf;

  p = inputParser;
  p.addOptional( 'outFile', defaultOutFile, @(x) true );
  p.addOptional( 'delayTime', defaultDelayTime, @isnumeric );
  p.addParameter('loopCount', defaultLoopCount, @isnumeric );
  p.parse( varargin{:} );
  outFile = p.Results.outFile;
  delayTime = p.Results.delayTime;
  loopCount = p.Results.loopCount;

  if ndims( cube ) ~= 3, error('cube must be a 3D array'); end

  first = 1;
  for i=1:size( cube, 3 )
    img = double( cube( :, :, i ) );
    img = uint8( img * 255 );

    if first == 1
      imwrite(img,outFile,'gif','LoopCount',loopCount,'DelayTime',delayTime);
      first = 0;
    else
      imwrite(img,outFile,'gif','WriteMode','append','DelayTime',delayTime);
    end
  end

end

