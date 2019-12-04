
function img = iwtDeaubechies2( wt, varargin )
  % img = iwtDeaubechies2( wt[, split] );
  % Performs an inverse Deaubechies wavelet transform of an image (with circular boundary conditions)
  %
  % Inputs:
  % img - 2D array representing the wavelet transform of the image
  %
  % Optional Inputs:
  % split - array specifying the number of levels of the wavelet transform.
  %   by default, split is 1 (indicating only one level).
  %   Example: [1 1; 1 0] will have 2 levels.  The size of the highest frequency
  %   portion in the final level will be double that of the other portions since
  %   it wasn't split.
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
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
  
  % H/L - High / Low pass filter
  % h/v - horizontal / vertical direction

  sWT = size(wt);
  wt11 = wt(1:sWT(1)/2,1:sWT(2)/2);
  wt21 = wt(sWT(1)/2+1:end,1:sWT(2)/2);
  wt12 = wt(1:sWT(1)/2,sWT(2)/2+1:end);
  wt22 = wt(sWT(1)/2+1:end,sWT(2)/2+1:end);

  sSplit = size( split );
  if max( mod( log2(sSplit), 1 ) ) ~= 0
    error( 'size of split should be even' );
  end

  nSplit = numel(split);
  if nSplit > 1
    sSplit = size(split);
    split11 = split(1:sSplit(1)/2,1:sSplit(2)/2);
    split12 = split(1:sSplit(1)/2,sSplit(2)/2+1:end);
    split21 = split(sSplit(2)/2+1:end,1:sSplit(1)/2);
    split22 = split(sSplit(2)/2+1:end,sSplit(2)/2+1:end);

    if sum(split11(:))>0
      wt11 = iwtDeaubechies2( wt11, split11 );
    end
    if sum(split12(:))>0
      wt12 = iwtDeaubechies2( wt12, split12 );
    end
    if sum(split21(:))>0
      wt21 = iwtDeaubechies2( wt21, split21 );
    end
    if sum(split22(:))>0
      wt22 = iwtDeaubechies2( wt22, split22 );
    end
  end


  tmp = upsample2( wt11, [1 2] );
  tmp3 = 3 * tmp;
  tmpSqrt3 = tmp * sqrt(3);
  wt1_1 = tmp - tmpSqrt3 + ...
    circshift( tmp3 - tmpSqrt3, [0 1] ) + ...
    circshift( tmp3 + tmpSqrt3, [0 2] ) + ...
    circshift( tmp + tmpSqrt3, [0 3] );

  tmp = upsample2( wt12, [1 2] );
  tmp3 = 3 * tmp;
  tmpSqrt3 = tmp * sqrt(3);
  wt1_2 = -( tmp + tmpSqrt3 ) + ...
    circshift( tmp3 + tmpSqrt3, [0 1] ) + ...
    circshift( -( tmp3 - tmpSqrt3 ), [0 2] ) + ...
    circshift( tmp - tmpSqrt3, [0 3] );

  wt1 = upsample2( wt1_1 + wt1_2, [2 1] );


  tmp = upsample2( wt21, [1 2] );
  tmp3 = 3 * tmp;
  tmpSqrt3 = tmp * sqrt(3);
  wt2_1 = tmp - tmpSqrt3 + ...
    circshift( tmp3 - tmpSqrt3, [0 1] ) + ...
    circshift( tmp3 + tmpSqrt3, [0 2] ) + ...
    circshift( tmp + tmpSqrt3, [0 3] );

  tmp = upsample2( wt22, [1 2] );
  tmp3 = 3 * tmp;
  tmpSqrt3 = tmp * sqrt(3);
  wt2_2 = -( tmp + tmpSqrt3 ) + ...
    circshift( tmp3 + tmpSqrt3, [0 1] ) + ...
    circshift( -( tmp3 - tmpSqrt3 ), [0 2] ) + ...
    circshift( tmp - tmpSqrt3, [0 3] );

  wt2 = upsample2( wt2_1 + wt2_2, [2 1] );


  tmp = wt1;
  tmp3 = 3 * tmp;
  tmpSqrt3 = tmp * sqrt(3);
  sig1 = tmp - tmpSqrt3 + ...
    circshift( tmp3 - tmpSqrt3, [1 0] ) + ...
    circshift( tmp3 + tmpSqrt3, [2 0] ) + ...
    circshift( tmp + tmpSqrt3, [3 0] );

  tmp = wt2;
  tmp3 = 3 * tmp;
  tmpSqrt3 = tmp * sqrt(3);
  sig2 = -( tmp + tmpSqrt3 ) + ...
    circshift( tmp3 + tmpSqrt3, [1 0] ) + ...
    circshift( -( tmp3 - tmpSqrt3 ), [2 0] ) + ...
    circshift( tmp - tmpSqrt3, [3 0] );
  
  img = ( sig1 + sig2 ) / 32;
end

