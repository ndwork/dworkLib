
function recon = reconNonRectSupports( supports, kDataEvenCols, kDataOddCols )
  % recon = reconNonRectSupport( supports, kDataEvenCols, kDataOddCols )
  %
  % Inputs:
  % supports - a two-dimensional binary mask with 1 indicating which pixels are in the support
  %            or a 3D array, one support for each coil
  % kDataEvenCols - 
  % kDataOddCols - 
  %
  % Outputs:
  % recon - a two dimensional image
  %
  % Written by Nicholas Dwork - Copyright 2022
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  nRows = size( supports, 1 );
  nCols = size( supports, 2 );
  nCoils = size( kDataEvenCols, 3 );

  if size( supports, 3 ) < nCoils
    supports = repmat( supports, [ 1 1 nCoils ] );
  end

  sFull = size( kDataEvenCols );
  sFull(2) = sFull(2) * 2;
  outerRows = outerRowsFromSupports( supports );
  innerRows = 1 - outerRows;
  nInnerRows = max( sum( innerRows, 1 ) );

  if size( outerRows, 2 ) > 1
    % Need shrink the outer rows of each column so that it has nInnerRows inner rows
    for coilIndx = 1 : nCoils
      coilInnerRows = 1 - outerRows(:,coilIndx);
      nCoilInnerRows = sum( coilInnerRows );
      if nCoilInnerRows < nInnerRows
        lastInnerRowIndx = find( coilInnerRows == 1, 1, 'last' );
        if numel( lastInnerRowIndx ) == 0, lastInnerRowIndx = numel( coilInnerRows ); end
        minIndx = max( 1, lastInnerRowIndx - nInnerRows + 1 );
        coilInnerRows( minIndx : lastInnerRowIndx ) = 1;
        if sum( coilInnerRows ) < nInnerRows
          coilInnerRows( 1 : nInnerRows ) = 1;
        end
        outerRows( :, coilIndx ) = 1 - coilInnerRows;
        innerRows( :, coilIndx ) = coilInnerRows;
      end
    end
  end

  kDataEvenZF = zeros( sFull );  % ZF - zero filled
  kDataEvenZF(:,2:2:end,:) = kDataEvenCols;
  coilReconsEvenCols = fftshift2( ifft2( ifftshift2( kDataEvenZF ) ) );

  outerSupports = zeros( size( supports ) );
  for coilIndx = 1 : nCoils
    outerSupports(:,:,coilIndx) = bsxfun( @times, supports(:,:,coilIndx), outerRows(:,coilIndx) );
  end

  if ismatrix( supports )
    coilReconsOuter = 2 * bsxfun( @times, coilReconsEvenCols, outerSupports );
  else
    coilReconsOuter = 2 * coilReconsEvenCols .* outerSupports;
  end

  kxFull = size2fftCoordinates( nCols );
  kxInnerMissing = kxFull(2:2:end);
  kyInner = size2fftCoordinates( nInnerRows );
  [ kxInnerMissingMesh, kyInnerMissingMesh ] = meshgrid( kxInnerMissing, kyInner );
  trajMissing = [ kyInnerMissingMesh(:) kxInnerMissingMesh(:) ];
  nxInnerMissing = numel( kxInnerMissing );
  kMissing = iGrid_2D( coilReconsEvenCols, trajMissing );
  kMissing = reshape( kMissing, [ nInnerRows nxInnerMissing nCoils ] );

  kInner = zeros( nInnerRows, nCols, nCoils );
  kInner(:,1:2:end,:) = kDataOddCols;
  kInner(:,2:2:end,:) = kMissing;

  % Now that we have kInner, we can subtract away the outer portion

  [ kxInnerMesh, kyInnerMesh ] = meshgrid( kxFull, kyInner );
  trajInner = [ kyInnerMesh(:) kxInnerMesh(:) ];

  tmp = iGrid_2D( coilReconsOuter, trajInner );
  kRemaining = kInner - reshape( tmp, [ nInnerRows nCols nCoils ] );

  coilReconsInner = grid_2D( reshape( kRemaining, [], nCoils ), trajInner, [ nRows nCols ] );
  for coilIndx = 1 : nCoils
    coilReconsInner(:,:,coilIndx) = bsxfun( @times, coilReconsInner(:,:,coilIndx), innerRows(:,coilIndx) );
  end

  coilRecons = coilReconsInner + coilReconsOuter;

  if nCoils > 1
    recon = mri_reconRoemer( coilRecons );
  else
    recon = coilRecons;
  end
end
