
function wt = wtHaar2( img, varargin )
  % wt = wtHaar2( sig [, split ] );
  % Performs a Haar wavelet transform of an image
  %
  % Inputs:
  % img - 2D array representing the image to be transformed
  %
  % Optional Inputs:
  % split - array specifying the number of levels of the wavelet transform.
  %   by default, split is 1 (indicating only one level).
  %   Example: [1 1; 1 0] will have 2 levels.  The size of the highest frequency
  %   portion in the final level will be double that of the other portions since
  %   it wasn't split.
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  defaultSplit = 1;
  p = inputParser;
  p.addOptional( 'split', defaultSplit );
  p.parse( varargin{:} );
  split = p.Results.split;

  sSplit = size( split );
  if max( mod( log2(sSplit), 1 ) ) ~= 0
    error( 'size of split should be even' );
  end

  wt1 = img(1:2:end-1,:) + img(2:2:end,:);
  wt11 = wt1(:,1:2:end-1) + wt1(:,2:2:end);
  wt12 = wt1(:,1:2:end-1) - wt1(:,2:2:end);

  wt2 = img(1:2:end-1,:) - img(2:2:end,:);
  wt21 = wt2(:,1:2:end-1) + wt2(:,2:2:end);
  wt22 = wt2(:,1:2:end-1) - wt2(:,2:2:end);

  if numel( split ) > 1
    sSplit = size(split);
    s11 = split( 1:sSplit(1)/2, 1:sSplit(2)/2 );
    s12 = split( 1:sSplit(1)/2, sSplit(2)/2+1:end );
    s21 = split( sSplit(2)/2+1:end, 1:sSplit(1)/2 );
    s22 = split( sSplit(2)/2+1:end, sSplit(2)/2+1:end );

    if sum( s11(:) ) > 0
      if max( mod(size(wt11),2) ) > 0
        error('wtHaar2: improper dimensions of image');
      end
      wt11 = wtHaar2( wt11, s11 );
    end
    if sum( s12(:) ) > 0
      if max( mod(size(wt12),2) ) > 0
        error('wtHaar2: improper dimensions of image');
      end
      wt12 = wtHaar2( wt12, s12 );
    end
    if sum( s21(:) ) > 0
      if max( mod(size(wt21),2) ) > 0
        error('wtHaar2: improper dimensions of image');
      end
      wt21 = wtHaar2( wt21, s21 );
    end
    if sum( s22(:) ) > 0
      if max( mod(size(wt22),2) ) > 0
        error('wtHaar2: improper dimensions of image');
      end
      wt22 = wtHaar2( wt22, s22 );
    end
  end

  wt = 0.5 * [ wt11 wt12; wt21 wt22 ];
end
