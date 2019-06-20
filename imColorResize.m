
function out = imColorResize( img, newSize, varargin )
  % out = imColorResize( img, newSize, varargin )
  %
  % Inputs:
  %   img - original color image
  %   newSize - 1, 2, or 3 element array specifying new array size;
  %     If 1 element, then assumed to be the scale
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

  sImg = size(img);
  if numel( newSize ) == 1
    newChannelSize = sImg(1:2) * newSize;
  else
    newChannelSize = newSize(1:2);
  end

  if ismatrix( img )
    nChannels = 1;
  else
    nChannels = sImg(3);
  end

  if isempty( gcp( 'nocreate' ) )
    out = zeros( [newChannelSize nChannels] );
    for ch = 1 : nChannels
      out(:,:,ch) = imresize( img(:,:,ch), newChannelSize, varargin{:} );
    end

  else
    out = cell( 1, 1, nChannels );
    parfor ch = 1 : nChannels
      out{ch} = imresize( img(:,:,ch), newChannelSize, varargin{:} );   %#ok<PFBNS>
    end
    out = cell2mat( out );

  end
end

