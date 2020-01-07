
function imH = wavShow( wt, varargin )
  % imH = wavshow( wt [, scale, 'split', split] );
  %
  % shows the wavelet transform with individual scaling of each portion for
  % improved viewing
  %
  % Inputs:
  % img - 2D array representing the wavelet transform of the image
  %
  % Optional Inputs:
  % split - array specifying the number of levels of the wavelet transform.
  %   by default, split is 1 (indicating only one level).
  %   Example: [1 1 1 0] will have 3 levels.  The size of the last portion
  %   in the final level will be double that of the other portions since it
  %   wasn't split.
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp( 'Usage:  imH = wavshow( wt [, scale, ''split'', split] )' );
    return
  end

  defaultSplit = 1;
  p = inputParser;
  p.addOptional( 'scale', 1, @ispositive );
  p.addParameter( 'range', [] );
  p.addParameter( 'split', defaultSplit );
  p.parse( varargin{:} );
  range = p.Results.range;
  scale = p.Results.scale;
  split = p.Results.split;

  if max( abs( imag( wt(:) ) ) ) ~= 0
    error( 'Input to wavShow is complex' );
  end
  
  img2show = wavScale( wt, split );
  imH = imshowscale( img2show, scale, 'range', range );
end

