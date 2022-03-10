
function traj = mri_makeTrajPts( nDim, type, varargin )
  % traj = makeTrajPts( nDim, type, parameters );
  %
  % traj is an N x nDim array, where N is the number of trajectory points
  %
  % Inputs:
  % nDim - the number of dimensions
  % type - the type of trajectory to create
  %        'poissonDisc' - poisson disc sampling
  %        'propeller'
  %        'random' - (default)
  %        'radial' - evenly spaced radial spokes
  %        'rosette' - rosette flower
  %        'spinWarp'
  %
  % Parameters for trajectories:
  %   'poissonDisc': traj = makeTrajPts( nDim, 'poissonDisc', radius );
  %     radius is the radius of the disc (nominal distance between points)
  %   'propeller': traj = makeTrajPts( nDim, nReadout, nLines, dkLine, ...
  %     nAngles );
  %   'random': traj = makeTrajPts( nDim, 'random', nTraj );
  %     nTraj is the number of points in the trajectory
  %   'radial': traj = makeTrajPts( nDim, 'radial', nSpokes, nPtsPerSpoke );
  %   'rosette': traj = makeTrajPts( nDim, 'rosette', nTraj, dt, f1, f2 );
  %     nTraj is the number of points in the trajectory
  %     dt is the spacing between samples
  %     f1 is the amplitude's frequency
  %     f2 is the phase's frequency
  %
  % To view the sampling pattern:
  %   plot( kTraj(:,1), kTraj(:,2), 'o', 'MarkerFaceColor', 'k', ...
  %     'MarkerEdgeColor', 'k', 'MarkerSize', 4 );
  %   set( gca, 'xTick', [], 'yTick', [] );
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin<2, type=''; end

  switch type
    case 'poissonDisc'
      traj = makeTrajPts_poissonDisc( nDim, varargin{:} );
    case 'propeller'
      traj = makeTrajPts_propeller( nDim, varargin{:} );
    case 'spinWarp'
      traj = makeTrajPts_spinWarp( nDim, varargin{:} );
    case 'random'
      traj = makeTrajPts_random(nDim, varargin{:} );
    case 'radial'
      traj = makeTrajPts_radial(nDim, varargin{:} );
    case 'rosette'
      traj = makeTrajPts_rosette(nDim, varargin{:} );
    otherwise
      traj = makeTrajPts_random(nDim, varargin{:} );
  end

  %plot( traj(:,1), traj(:,2), 'o', 'MarkerFaceColor', 'k', ...
  %  'MarkerEdgeColor', 'k', 'MarkerSize', 2 );
  %set( gca, 'xTick', [], 'yTick', [] );
end


function traj = makeTrajPts_poissonDisc( nDim, gamma )
  if nDim ~= 2, error('Not yet implemented'); end
  
  Delta = 0.15;
  r = @(x) ( norms( x, 2, 2 ) + Delta ) / gamma;
  min_r = r( [0 0] );
  traj = poissonDisc2( r, 'min_r', min_r, 'bounds', bounds );
end


function traj = makeTrajPts_propeller( nDim, nReadout, nLines, dkLine, nAngles )
  if nDim ~= 2, error('Propeller is a 2D Trajectory'); end

  nKperAngle = nReadout * nLines;
  dAngle = pi / nAngles;

  kx = ones(nLines,1) * linspace(-0.5,0.5,nReadout);
  kyExtent = dkLine * (nLines-1);
  kyMax = kyExtent/2;
  ky = (-kyMax:dkLine:kyMax)' * ones(1,nReadout);

  traj = zeros(nReadout*nLines*nAngles,2);
  for i=1:nAngles
    thisAngle = (i-1)*dAngle;
    thisKx = cos(thisAngle)*kx(:) - sin(thisAngle)*ky(:);
    thisKy = sin(thisAngle)*kx(:) + cos(thisAngle)*ky(:);

    traj((i-1)*nKperAngle+1:i*nKperAngle,1) = thisKx;
    traj((i-1)*nKperAngle+1:i*nKperAngle,2) = thisKy;
  end

  trajNorms = norms( traj, 2, 2 );
  traj = traj * 0.5 / max( abs( trajNorms ) );
end


function traj = makeTrajPts_spinWarp( nDim, dk )
  k = transpose( -0.5 : dk : 0.5-dk );
  traj = zeros( numel(k)^nDim, nDim );

  if nDim == 1
    traj = k;
  elseif nDim == 2
    [kx,ky] = meshgrid( k, k );
    traj(:,1) = kx(:);
    traj(:,2) = ky(:);
  elseif nDim==3
    [kx,ky,kz] = meshgrid( k, k, k );
    traj(:,1) = kx(:);
    traj(:,2) = ky(:);
    traj(:,3) = kz(:);
  end
end


function traj = makeTrajPts_random( nDim, nTraj )
  traj = rand( nTraj, nDim ) - 0.5;
end


function traj = makeTrajPts_radial( nDim, nSpokes, nPtsPerSpoke )
  if nDim ~= 2, error('Radial trajectory just for 2 dimensions'); end

  thetas = linspace( 0, 2*pi, nSpokes + 1 );
  thetas = thetas( 1 : end-1 );
  rs = linspace( 0, 0.5, nPtsPerSpoke );

  traj = zeros( nSpokes*nPtsPerSpoke, 2 );
  for i=1:nSpokes
    lowIndx = (i-1)*nPtsPerSpoke+1;
    highIndx = i*nPtsPerSpoke;
    traj(lowIndx:highIndx,1) = rs * cos(thetas(i));
    traj(lowIndx:highIndx,2) = rs * sin(thetas(i));
  end
end


function traj = makeTrajPts_rosette( nDim, nPtsPerCycle, nCycles, f1, f2 )
  % Made according to "Multishot Rosette Trajectories for Spectrally Selective MR Imaging" by Noll., 1997
  if nDim ~= 2, error('Rosette trajectory just for 2 dimensions'); end

  nPts = nPtsPerCycle * nCycles;
  cyclePts = zeros( nPtsPerCycle, 2 );
  traj = zeros( nPts, 2 );

  t = ( 0 : nPtsPerCycle-1 ) / nPtsPerCycle;

  eIndx = 0;
  for cycle = 1 : nCycles
    dOffset = ( pi / nCycles ) * ( cycle - 1 );

    cyclePts(:,1) = 0.5 * cos( 2*pi * f1 * t + dOffset ) .* sin( 2*pi * f2 * t + dOffset );
    cyclePts(:,2) = 0.5 * cos( 2*pi * f1 * t + dOffset ) .* cos( 2*pi * f2 * t + dOffset );

    sIndx = eIndx + 1;
    eIndx = sIndx + nPtsPerCycle - 1;
    traj( sIndx : eIndx, : ) = cyclePts;
  end

end



