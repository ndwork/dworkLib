
function out = makeAnaglyph( imgL, imgR )
  % This function makes an anaglyph from two individual images
  %
  % out = makeAnaglyph( imgL, imgR )
  %
  % Inputs:
  % imgL - the image that will become red
  % imgR - the image that will become blue
  %
  % Outputs:
  % out - the anaglyph image
  %
  % Written by Nicholas Dwork - Copyright 2017
  %
  % This software is offered under the GNU General Public License 3.0.  It 
  % is offered without any warranty expressed or implied, including the 
  % implied warranties of merchantability or fitness for a particular 
  % purpose.

  if ~ismatrix( imgL ), imgL = rgb2gray(imgL); end;
  if ~ismatrix( imgR ), imgR = rgb2gray(imgR); end;

  sImgL = size( imgL );
  out = zeros( sImgL(1), sImgL(2), 3 );
  
  out(:,:,1) = imgL;
  out(:,:,3) = imgR;
end