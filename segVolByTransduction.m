
function out = segVolByTransduction( vol, ins, outs, varargin )
  % out = segVolByTransduction( vol, ins, outs [, 'lambda', lambda, 'nbhdSize', nbhdSize ] )
  %
  % Implements the algorithm of "Segmentation by Transduction" by Duchenne et al. for a
  % three-dimensional volume
  %
  % Inputs:
  % vol - a three-dimensional array
  % ins - a 1D array of indices into img for pixels that are within the segmented region
  % outs - a 1D array of indices into img for pixels that are outside of the segmented region
  %
  % Optional Inputs:
  % nbhdSize - the neighborhood size for those pixels that are considered local
  %
  % Outputs:
  % out - an array of size(vol) where pixels within/outside of the segmentation have value 1/0
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

  if numel( nbhdSize ) == 1, nbhdSize = [ nbhdSize nbhdSize nbhdSize ]; end
  if numel( nbhdSize ) == 2, error( 'Neighborhood is a three element value' ); end
  sVol = size( vol );
  nVol = numel( vol );

  knownMask = zeros( size( vol, 1 ), size( vol, 2 ), size( vol, 3 ) );
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
  nbhdYs = nbhdCoords{1};
  nbhdXs = nbhdCoords{2};
  nbhdZs = nbhdCoords{3};
  [ nbhdYsTmp, nbhdXsTmp, nbhdZsTmp ] = ndgrid( nbhdYs, nbhdXs, nbhdZs );
  nbhdDists = sqrt( nbhdYsTmp.^2 + nbhdXsTmp.^2 + nbhdZsTmp.^2 );
  dWeights = exp( -lambda * nbhdDists );  %distance Weights
  clear nbhdYsTmp nbhdYsTmp nbhdZsTmp

  aIndx = 1;
  aRowIndx = 1;
  for k = 1 : sVol(3)
    disp(['Working on k=', indx2str(k,sVol(3)), ' of ', num2str(sVol(3)) ]);
    for i = 1 : sVol(2)
      for j = 1 : sVol(1)
  
        if knownMask(j,i,k) ~= 0
          aRows(aIndx) = aRowIndx;
          aVals(aIndx) = 1;
          aCols(aIndx) = sub2ind(sVol,j,i,k);
          b(aRowIndx) = knownMask(j,i,k);
          aRowIndx = aRowIndx + 1;
          aIndx = aIndx + 1;
          continue;
        end
  
        for iNbhd = nbhdXs(:)'
          if i + iNbhd < 1 || i + iNbhd > sVol(2), continue; end
          iIndxD = iNbhd - min(nbhdXs(:)) + 1;
  
          for jNbhd = nbhdYs(:)'
            if j + jNbhd < 1 || j + jNbhd > sVol(1), continue; end
            jIndxD = jNbhd - min(nbhdYs(:)) + 1;
  
            for kNbhd = nbhdZs(:)'
              if k + kNbhd < 1 || k + kNbhd > sVol(3), continue; end
              if kNbhd == 0 && jNbhd == 0 && iNbhd == 0, continue; end
              kIndxD = kNbhd - min(nbhdZs(:)) + 1;
  
  
              cImg = vol(j,i,k);
              cNbhd = vol(j+jNbhd,i+iNbhd,k+kNbhd);
              w = dWeights(jIndxD,iIndxD,kIndxD) * exp( -norm( cImg(:) - cNbhd(:) ) );
    
              if knownMask(j+jNbhd,i+iNbhd,k+kNbhd) ~= 0
                aRows(aIndx) = aRowIndx;  aCols(aIndx) = sub2ind(sVol,j,i,k);
                aVals(aIndx) = w;
                b(aRowIndx) = w * knownMask(j+jNbhd,i+iNbhd,k+kNbhd);
                aIndx = aIndx + 1;
              else
                aRows( aIndx ) = aRowIndx;
                aVals( aIndx ) =  w;
                aCols( aIndx ) = sub2ind( sVol, j, i, k );  
                aIndx = aIndx + 1;
                aRows( aIndx ) = aRowIndx;
                aVals( aIndx ) = -w;
                aCols( aIndx ) = sub2ind( sVol, j+jNbhd, i+iNbhd, k+kNbhd );
                aIndx = aIndx + 1;
              end

              aRowIndx = aRowIndx + 1;
            end
          end
        end
      end
    end
  end
  aRows = aRows( 1 : aIndx-1 );
  aCols = aCols( 1 : aIndx-1 );
  aVals = aVals( 1 : aIndx-1 );
  b = b( 1 : max(aRows) );

  A = sparse( aRows, aCols, aVals, max(aRows), nVol );

  out = A \ b;
  out = reshape( out > 0, sVol );
end





