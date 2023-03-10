
function out = segImgByTransduction( img, ins, outs, varargin )
  % out = segByTransduction( img, ins, outs [, 'lambda', lambda, 'nbhdSize', nbhdSize ] )
  %
  % Implements the algorithm of "Segmentation by Transduction" by Duchenne et al.
  %
  % Inputs:
  % img - an array
  % ins - a 1D array of indices into img for pixels that are within the segmented region
  % outs - a 1D array of indices into img for pixels that are outside of the segmented region
  %
  % Optional Inputs:
  % nbhdSize - the neighborhood size for those pixels that are considered local
  %
  % Outputs:
  % out - an array of [ nRows x nCols ] where pixels within/outside of the
  %       segmentation have value 1/0
  %
  % Written by Nicholas Dwork - Copyright 2022
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addParameter( 'lambda', 1, @isnonnegative );
  p.addParameter( 'nbhdSize', 5, @ispositive );
  p.parse( varargin{:} );
  lambda = p.Results.lambda;
  nbhdSize = p.Results.nbhdSize;

  if numel( nbhdSize ) == 1, nbhdSize = [ nbhdSize nbhdSize ]; end

  for dim = 1 : size( img, 3 )
    dimImg = img(:,:,dim);
    img(:,:,dim) = ( dimImg - mean( dimImg(:) ) ) / std( dimImg(:) );
  end
  if ismatrix(img), img = repmat( img, [ 1 1 3 ] ); end  % convert into gray color image

  sImg = size( img );
  nImg = sImg(1) * sImg(2);

  knownMask = zeros( size( img, 1 ), size( img, 2 ) );
  knownMask(ins) = 1;
  knownMask(outs) = -1;
  nKnown = sum( knownMask(:) ~= 0 );
  nNotKnown = sum( knownMask(:) == 0 );

  aRows = zeros( 2 * nNotKnown * ( prod( nbhdSize ) - 1 ) + nKnown, 1 );
  aCols = zeros( 2 * nNotKnown * ( prod( nbhdSize ) - 1 ) + nKnown, 1 );
  aVals = zeros( 2 * nNotKnown * ( prod( nbhdSize ) - 1 ) + nKnown, 1 );
  b     = zeros( nNotKnown * ( prod( nbhdSize ) - 1 ) + nKnown, 1 );

  % Find the spatial distance weights for those points in the local neighborhood
  nbhdCoords = size2imgCoordinates( nbhdSize );
  nbhdYs = nbhdCoords{1};  sdevY = std( size2imgCoordinates( size(img,1) ) );
  nbhdXs = nbhdCoords{2};  sdevX = std( size2imgCoordinates( size(img,2) ) );
  [ nbhdYsTmp, nbhdXsTmp ] = ndgrid( nbhdYs / sdevY, nbhdXs / sdevX );
  nbhdDistsSq = nbhdYsTmp .* nbhdYsTmp + nbhdXsTmp .* nbhdXsTmp;
  clear nbhdYsTmp nbhdYsTmp

  aIndx = 1;
  aRowIndx = 1;
  for i = 1 : sImg(2)
    for j = 1 : sImg(1)

      if knownMask(j,i) ~= 0
        aRows(aIndx) = aRowIndx;  aVals(aIndx) = 1;  aCols(aIndx) = sub2ind(sImg(1:2),j,i);
        b(aRowIndx) = knownMask(j,i);
        aRowIndx = aRowIndx + 1;
        aIndx = aIndx + 1;
        continue;
      end

      cImg = img(j,i,:);

      for iNbhd = nbhdXs(:)'
        if i + iNbhd < 1 || i + iNbhd > sImg(2), continue; end
        iIndxD = iNbhd - min(nbhdXs(:)) + 1;

        for jNbhd = nbhdYs(:)'
          if j + jNbhd < 1 || j + jNbhd > sImg(1), continue; end
          if iNbhd == 0 && jNbhd == 0, continue; end
          jIndxD = jNbhd - min(nbhdYs(:)) + 1;

          nbhdDistSq = nbhdDistsSq( jIndxD, iIndxD );
          cNbhd = img( j+jNbhd, i+iNbhd, : );
          cDistSq = norm( cImg(:) - cNbhd(:) )^2;
          w = exp( -lambda * sqrt( nbhdDistSq + cDistSq ) );

          if knownMask(j+jNbhd,i+iNbhd) ~= 0
            aRows(aIndx) = aRowIndx;  aCols(aIndx) = sub2ind(sImg(1:2),j,i);
            aVals(aIndx) = w;
            b(aRowIndx) = w * knownMask(j+jNbhd,i+iNbhd);
            aIndx = aIndx + 1;
          else
            aRows(aIndx) = aRowIndx;  aVals(aIndx) =  w;  aCols(aIndx) = sub2ind(sImg(1:2),j,i);  
            aIndx = aIndx + 1;
            aRows(aIndx) = aRowIndx;  aVals(aIndx) = -w;  aCols(aIndx) = sub2ind(sImg(1:2),j+jNbhd,i+iNbhd);
            aIndx = aIndx + 1;
          end

          aRowIndx = aRowIndx + 1;
        end
      end
    end
  end
  aRows = aRows( 1 : aIndx-1 );
  aCols = aCols( 1 : aIndx-1 );
  aVals = aVals( 1 : aIndx-1 );
  b = b( 1 : max(aRows) );

  A = sparse( aRows, aCols, aVals, max(aRows), nImg );

  out = A \ b;
  out = reshape( out > 0, sImg(1:2) );
end





