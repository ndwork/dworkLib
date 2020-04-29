
function [griddedPts,samples] = movePointsToGrid( pts, mins, maxs, Ns )
  % [griddedPts,samples] = movePointsToGrid( pts, mins, maxs, nums )
  %
  % Moves points in an Real^N space to the points specified by grid
  %
  % Inputs:
  % pts - a NxM array (where M is the number of points)
  % mins - a 1D array of size N specifying the minimum coordinate of each dimension of the grid
  % maxs - a 1D array of size N specifying the maximum coordinate of each dimension of the grid
  % Ns - a 1D array of size N specifying the number of grid points along each dimension
  %
  % Outputs:
  % griddedPts - an NxM array where the point coordinates have been shifted so that they lie on the grid
  %
  % Optional Outputs:
  % samples - a array of size Ns where each element corresponds to a grid point and the value of 
  %   each bin is the number of samples on that grid point
  %
  % Written by Nicholas Dwork, Copyright 2019
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  griddedPts = zeros( size(pts) );

  nDims = numel( mins );
  if numel( maxs ) ~= nDims, error('maxs must have the same size as mins'); end
  if numel( Ns ) ~= nDims, error('Ns must have the same size as mins'); end

  griddedCoords = cell( nDims, 1 );
  for j=1:nDims
    griddedCoords{j} = linspace( mins(j), maxs(j), Ns(j) );
  end

  if nargout > 1, samples = zeros( Ns(:)' ); end  % Create a sample array

  for i = 1:size( pts, 2 )

    minIndxs = cell(nDims,1);
    for j=1:nDims
      %disp([ '(j,i): ', num2str(j), ', ', num2str(i) ]);

      theseCoords = griddedCoords{j};
      [~,minIndx] = min( abs( pts(j,i) - theseCoords ) );
      griddedPts(j,i) = theseCoords( minIndx );

      minIndxs{j} = minIndx;
    end

    if nargout > 1
      samples( minIndxs{:} ) = samples( minIndxs{:} ) + 1;
    end
  end

end
