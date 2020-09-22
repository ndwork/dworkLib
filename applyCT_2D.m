

function out = applyCT_2D( fftData, traj, N, kCy, kCx, Cy, Cx, varargin )
  % out = applyCT_2D( fftData, traj, N, kCy, kCx, Cy, Cx [, gridKs, 'type', type ] )
  %
  % Inputs:
  % fftData - Mx1 array specifying Fourier value at traj(indx,:)
  % traj - Mx2 array where first/second column represents ky/kx location
  %   of trajectory
  % N - 2 element array specifying [Ny Nx] points in grid
  %
  % Written by Nicholas Dwork - Copyright 2016

  defaultGridKs = [];
  defaultType = [];
  p = inputParser;
  p.addOptional( 'gridKs', defaultGridKs );
  p.addParameter( 'type', defaultType );
  p.parse( varargin{:} );
  gridKs = p.Results.gridKs;
  type = p.Results.type;

  if numel( gridKs ) == 0
    gridKs = size2fftCoordinates( N );
    gridKy=gridKs{1};  gridKx=gridKs{2};
    [gridKx,gridKy] = meshgrid(gridKx,gridKy);
  else
    gridKy = gridKs(:,1);
    gridKx = gridKs(:,2);
  end

  kws = [ max(kCy), max(kCx) ];
  kDistThreshY = kws(1);
  kDistThreshX = kws(2);

  nTraj = size( traj, 1 );
  out = zeros( nTraj, 1 );
  parfor trajIndx = 1:nTraj
    kDistsY = abs( traj(trajIndx,1) - gridKy );
    kDistsX = abs( traj(trajIndx,2) - gridKx );  %#ok<PFBNS>
    shortDistIndxs = find( kDistsY < kDistThreshY & kDistsX < kDistThreshX );
    shortDistsKy = kDistsY( shortDistIndxs );
    shortDistsKx = kDistsX( shortDistIndxs );
    CValsY = interp1( kCy, Cy, shortDistsKy, 'linear', 0 );
    CValsX = interp1( kCx, Cx, shortDistsKx, 'linear', 0 );
    kVals = fftData( shortDistIndxs );                                           %#ok<PFBNS>
    out( trajIndx ) = sum( kVals .* CValsY .* CValsX );
  end

  onesCol = ones(nTraj,1);
  if ~strcmp( type, 'noCirc' )
    for dim=1:2
      alt = zeros( size(traj) );
      alt(:,dim) = onesCol;

      for altDir=[-1 1]
        NewTraj = traj + altDir*alt;
        if altDir < 0
          NewTrajIndxs = find( NewTraj(:,dim) > -0.5-kws(dim) );
        else
          NewTrajIndxs = find( NewTraj(:,dim) < 0.5+kws(dim) );
        end

        NewTraj = NewTraj( NewTrajIndxs, : );
        for i=1:numel(NewTrajIndxs)
          trajIndx = NewTrajIndxs(i);
          NewDistsKy = abs( NewTraj(i,1) - gridKy );
          NewDistsKx = abs( NewTraj(i,2) - gridKx );
          NewShortDistIndxs = find( NewDistsKy < kDistThreshY & ...
                                    NewDistsKx < kDistThreshX );
          NewShortDistsKy = NewDistsKy( NewShortDistIndxs );
          NewShortDistsKx = NewDistsKx( NewShortDistIndxs );
          NewCValsY = interp1( kCy, Cy, NewShortDistsKy, 'linear', 0 );
          NewCValsX = interp1( kCx, Cx, NewShortDistsKx, 'linear', 0 );
          NewKVals = fftData( NewShortDistIndxs );
          out(trajIndx) = out(trajIndx) + sum( NewKVals .* NewCValsY .* NewCValsX );
        end
      end
    end
  end

end

