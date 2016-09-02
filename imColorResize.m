
function out = imColorResize( img, newSize, varargin )
  % out = imColorResize( img, newSize, varargin )
  %
  % Inputs:
  %   img - original color image
  %   newSize - 3 element array specifying new array size;
  %     3rd element must be 3.  If two elements are supplied, assumed 3rd
  %     element is 3.
  % 
  % Optional Inputs:
  %   All optional inputs accepted by imresize are accepted
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  newChannelSize = newSize(1:2);
  
  out = zeros( [newChannelSize 3] );
  
  out(:,:,1) = imresize( img(:,:,1), newChannelSize, varargin{:} );
  out(:,:,2) = imresize( img(:,:,2), newChannelSize, varargin{:} );
  out(:,:,3) = imresize( img(:,:,3), newChannelSize, varargin{:} );
end

