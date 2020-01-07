
function wt = wtDeaubechies2( img, varargin )
  % wt = wtDeaubechies2( sig [, split] );
  %
  % Performs a Deaubechies wavelet transform of an image (with circular boundary conditions)
  % Based on the Wikipedia page on the Daubechies wavelet transform and
  % (http://wavelets.pybytes.com/wavelet/db2/)
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

  sSplit = size( split );
  if max( mod( log2(sSplit), 1 ) ) ~= 0
    error( 'size of split should be even' );
  end

  imgSqrt3 = img * sqrt(3);
  img3 = 3 * img;

  imgPimgSqrt3 = img + imgSqrt3;
  img3PimgSqrt3 = img3 + imgSqrt3;
  imgMimgSqrt3 = img - imgSqrt3;
  img3MimgSqrt3 = img3 - imgSqrt3;
  
  wt1 = ( imgMimgSqrt3 + ...
    circshift( img3MimgSqrt3, [ -1 0 ] ) + ...
    circshift( img3PimgSqrt3, [ -2 0 ] ) + ...
    circshift( imgPimgSqrt3, [ -3 0 ] ) );
  wt1 = wt1(1:2:end,:);

  wt1Sqrt3 = wt1 * sqrt(3);
  wt13 = 3 * wt1;

  wt1Pwt1Sqrt3 = wt1 + wt1Sqrt3;
  wt13Pwt1Sqrt3 = wt13 + wt1Sqrt3;
  wt1Mwt1Sqrt3 = wt1 - wt1Sqrt3;
  wt13Mwt1Sqrt3 = wt13 - wt1Sqrt3;

  wt11 = ( wt1Mwt1Sqrt3 + ...
    circshift( wt13Mwt1Sqrt3, [ 0 -1 ] ) + ...
    circshift( wt13Pwt1Sqrt3, [ 0 -2 ] ) + ...
    circshift( wt1Pwt1Sqrt3, [ 0 -3 ] ) );
  wt11 = wt11(:,1:2:end);

  wt12 = ( -wt1Pwt1Sqrt3 + ...
    circshift( wt13Pwt1Sqrt3, [ 0 -1 ] ) + ...
    circshift( -wt13Mwt1Sqrt3, [ 0 -2 ] ) + ...
    circshift( wt1Mwt1Sqrt3, [ 0 -3 ] ) );
  wt12 = wt12(:,1:2:end);

  wt2 = ( -imgPimgSqrt3 + ...
    circshift( img3PimgSqrt3, [ -1 0 ] ) + ...
    circshift( -img3MimgSqrt3, [ -2 0 ] ) + ...
    circshift( imgMimgSqrt3, [ -3 0 ] ) );
  wt2 = wt2(1:2:end,:);

  wt2Sqrt3 = wt2 * sqrt(3);
  wt23 = 3 * wt2;

  wt2Pwt2Sqrt3 = wt2 + wt2Sqrt3;
  wt23Pwt2Sqrt3 = wt23 + wt2Sqrt3;
  wt2Mwt2Sqrt3 = wt2 - wt2Sqrt3;
  wt23Mwt2Sqrt3 = wt23 - wt2Sqrt3;

  wt21 = ( wt2Mwt2Sqrt3 + ...
    circshift( wt23Mwt2Sqrt3, [ 0 -1 ] ) + ...
    circshift( wt23Pwt2Sqrt3, [ 0 -2 ] ) + ...
    circshift( wt2Pwt2Sqrt3, [ 0 -3 ] ) );
  wt21 = wt21(:,1:2:end);

  wt22 = ( -wt2Pwt2Sqrt3 + ...
    circshift( wt23Pwt2Sqrt3, [ 0 -1 ] ) + ...
    circshift( -wt23Mwt2Sqrt3, [ 0 -2 ] ) + ...
    circshift( wt2Mwt2Sqrt3, [ 0 -3 ] ) );
  wt22 = wt22(:,1:2:end);


  nSplit = numel(split);
  if nSplit > 1
    sSplit = size(split);
    s11 = split(1:sSplit(1)/2,1:sSplit(2)/2);
    s12 = split(1:sSplit(1)/2,sSplit(2)/2+1:end);
    s21 = split(sSplit(2)/2+1:end,1:sSplit(1)/2);
    s22 = split(sSplit(2)/2+1:end,sSplit(2)/2+1:end);

    if sum(s11(:))>0
      if max( mod(size(wt11),2) ) > 0
        error('wtDeaubechies2: improper dimensions of image');
      end
      wt11 = wtDeaubechies2( wt11, s11 );
    end
    if sum(s12(:))>0
      if max( mod(size(wt12),2) ) > 0
        error('wtDeaubechies2: improper dimensions of image');
      end
      wt12 = wtDeaubechies2( wt12, s12 );
    end
    if sum(s21(:))>0
      if max( mod(size(wt21),2) ) > 0
        error('wtDeaubechies2: improper dimensions of image');
      end
      wt21 = wtDeaubechies2( wt21, s21 );
    end
    if sum(s22(:))>0
      if max( mod(size(wt22),2) ) > 0
        error('wtDeaubechies2: improper dimensions of image');
      end
      wt22 = wtDeaubechies2( wt22, s22 );
    end
  end

  wt = [ wt11 wt12; wt21 wt22 ] / 32;
end



