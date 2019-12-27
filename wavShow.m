
function imH = wavShow( wt, varargin )
  % imH = wavshow( sig [, scale, 'split', split] );
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


function out = wavScale( wt, split )

  sWT = size( wt );
  wt11 = wt( 1 : sWT(1)/2, 1 : sWT(2)/2 );
  wt12 = wt( 1 : sWT(1)/2, sWT(2)/2+1 : end );
  wt21 = wt( sWT(1)/2+1 : end, 1 : sWT(2)/2 );
  wt22 = wt( sWT(1)/2+1 : end, sWT(2)/2+1 : end );

  nSplit = numel(split);
  if nSplit > 1
    sSplit = size(split);
    s11 = split( 1:sSplit(1)/2, 1:sSplit(2)/2 );
    s12 = split( 1:sSplit(1)/2, sSplit(2)/2+1:end );
    s21 = split( sSplit(2)/2+1:end, 1:sSplit(1)/2 );
    s22 = split( sSplit(2)/2+1:end, sSplit(2)/2+1:end );

    if sum( s11(:) ) > 0
      if max( mod(size(wt11),2) ) > 0
        error('wavShow: improper dimensions of image');
      end
      wt11 = wavScale( wt11, s11 );
    end
    if sum( s12(:) ) > 0
      if max( mod(size(wt12),2) ) > 0
        error('wavShow: improper dimensions of image');
      end
      wt12 = wavScale( wt12, s12 );
    end
    if sum( s21(:) ) > 0
      if max( mod(size(wt21),2) ) > 0
        error('wavShow: improper dimensions of image');
      end
      wt21 = wavScale( wt21, s21 );
    end
    if sum( s22(:) ) > 0
      if max( mod(size(wt22),2) ) > 0
        error('wavShow: improper dimensions of image');
      end
      wt22 = wavScale( wt22, s22 );
    end
  end

  scaled11 = scaleImg( wt11, [0 1] );
  scaled12 = scaleImg( wt12, [0 1] );
  scaled21 = scaleImg( wt21, [0 1] );
  scaled22 = scaleImg( wt22, [0 1] );

  out = [ scaled11 scaled12; scaled21 scaled22 ];
end

