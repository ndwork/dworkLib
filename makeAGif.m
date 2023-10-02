
function makeAGif( varargin )
  % makeAGif( [ inDir, outFile, delayTime, 'loopCount', loopCount ] )
  %
  %   This function makes a gif out of a set of images
  %   It assumes that all images are stored in their own directory (and
  %   no other files are stored in the directory with them).
  %   The images must be named so that they are listed in the order of display
  %   delayTime is the time between frames in the gif animation.
  %
  %   Inputs: 
  %   inDir (optional) - a string containing the name of the directory with input images
  %     By default, inDir is assumed to be './inDir'
  %   outFile (optional) - a string contianing the name of the gif file
  %     By default, outFiles is assumed ot be 'makeAGifOut.gif'
  %   delayTime (optional) - the time between frames in the animation
  %     By default, delayTime is assumed to be 0.3 seconds
  %   loopCount (optional) - the number of times to loop the animation
  %     By default, loopCount is infinity (meaning an endless loop)
  %
  % Written by Nicholas Dwork - Copyright 2015
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  defaultInDir = './inDir';
  defaultOutFile = 'makeAGifOut.gif';
  defaultDelayTime = 0.3;
  defaultLoopCount = Inf;

  p = inputParser;
  p.addOptional( 'inDir', defaultInDir, @(x) true );
  p.addOptional( 'outFile', defaultOutFile, @(x) true );
  p.addOptional( 'delayTime', defaultDelayTime, @isnumeric );
  p.addParameter('loopCount', defaultLoopCount, @isnumeric );
  p.parse( varargin{:} );
  inDir = p.Results.inDir;
  outFile = p.Results.outFile;
  delayTime = p.Results.delayTime;
  loopCount = p.Results.loopCount;

  if numel( outFile ) == 0, outFile = defaultOutFile; end

  imgs = dir( inDir );
  imgs = imgs(3:end);
  if numel(imgs)<1, disp('No images found'); end

  first = 1;
  for i=1:numel(imgs)
    if regexp( imgs(i).name, '^\.' ), continue; end
    if regexp( imgs(i).name, '\.gif$' ), continue; end

    img = imread( [inDir,'/',imgs(i).name] );
    if ismatrix(img)
      img = repmat(img,[1,1,3]);
    end
    [A,map] = rgb2ind(img,256);

    if first == 1
      imwrite(A,map,outFile,'gif','LoopCount',loopCount,'DelayTime',delayTime);
      first = 0;
    else
      imwrite(A,map,outFile,'gif','WriteMode','append','DelayTime',delayTime);
    end
  end

end
