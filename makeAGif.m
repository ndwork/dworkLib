
function makeAGif( varargin )
  % makeAGif( inDir, outFile, delayTime )
  %   This function makes a gif out of a set of images
  %   It assumes that all images are stored in their own directory (and
  %   no other files are stored in the directory with them).
  %   The images must be named so that they are listed in the order of display
  %   delayTime is the time between frames in the gif animation.
  %   inDir (optional) - a string containing the name of the directory with input images
  %     By default, inDir is assumed to be './inDir'
  %   outFile (optional) - a string contianing the name of the gif file
  %     By default, outFiles is assumed ot be 'makeAGifOut.gif'
  %   delayTime (optional) - the time between frames in the animation
  %     By default, delayTime is assumed to be 0.3 seconds
  %
  % Written by Nicholas Dwork - (c) 2015


  defaultInDir = './inDir';
  defaultOutFile = 'makeAGifOut.gif';
  defaultDelayTime = 0.3;

  p = inputParser;
  p.addOptional('inDir',defaultInDir);
  p.addOptional('outFile',defaultOutFile);
  p.addOptional('delayTime',defaultDelayTime,@isnumeric);
  p.parse( varargin{:} );
  inDir = p.Results.inDir;
  outFile = p.Results.outFile;
  delayTime = p.Results.delayTime;


  imgs = dir( inDir );
  imgs = imgs(3:end);
  if numel(imgs)<1, disp('No images found'); end;

  for i=1:numel(imgs)
    img = imread( [inDir,'/',imgs(i).name] );
    if ndims(img)==2
      img = repmat(img,[1,1,3]);
    end
    [A,map] = rgb2ind(img,256);

    if i == 1;
      imwrite(A,map,outFile,'gif','LoopCount',Inf,'DelayTime',delayTime);
    else
      imwrite(A,map,outFile,'gif','WriteMode','append','DelayTime',delayTime);
    end
  end

end
