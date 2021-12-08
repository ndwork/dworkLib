
function [ ks, c ] = findRadialDistortionFromLines( lines, varargin )
  % k = findRadialDistortionFromLines( img, lines [, order ] )
  %
  % Note: this function assumes the aspect ratio of point coordinates is 1.
  %
  % Inputs:
  % lines - a 1D cell array where each cell contains a set of points belonging to a line
  %   The element of each cell is a NL x 2 array where NL is the number of points for that line
  %
  % Optional Inputs:
  % order - the order of the parameters to identify
  %
  % Outputs:
  % ks - a 1D array specifying the radial distortion parameters
  % c - a 2 element array specifying the center of the image
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addOptional( 'order', 1, @ispositive );
  p.parse( varargin{:} );
  order = p.Results.order;

  ks0 = zeros( order, 1 );
  c0 = [0; 0;];

  %err = computeLineError( lines{1}, ks0, c0 );

  function err = computeThisError( x )
    %cThis = x(1:2);
    %ksThis = x(3:end);
    cThis = [0; 0;];
    ksThis = x;
    err = computeError( lines, ksThis, cThis );
  end

  %x0 = [ c0; ks0; ];
  x0 = ks0;
  [x,finalErr] = fminsearch( @computeThisError, x0 );   %#ok<ASGLU>

  %c = x(1:2);
  %ks = x(3:end);
  ks = x;
  c = [0; 0;];
end


function err = computeError( lines, ks, c )
  err = 0;
  for i = 1 : numel( lines )
    err = err + computeLineError( lines{i}, ks, c ).^2;
  end
  err = sqrt( err );
end


function err = computeLineError( pts, ks, c )
  rs = sqrt( pts(:,1).^2 + pts(:,2).^2 );
  rPower = rs;
  Ls = ones( size( rs ) );
  for kIndx = 1 : numel( ks )
    Ls = Ls + ks(kIndx) .* rPower;
    if kIndx < numel( ks ), rPower = rPower .* rs; end
  end

  undistortedPts = pts;
  undistortedPts(:,1) = c(1) + ( pts(:,1) - c(1) ) .* Ls;
  undistortedPts(:,2) = c(2) + ( pts(:,2) - c(2) ) .* Ls;

  % normalize the points
  %undistortedPts = normalizePts2D( undistortedPts );

  %figure; plotnice( pts(:,1), pts(:,2) );
  %hold all;  plotnice( undistortedPts(:,1), undistortedPts(:,2) )

  line.pt = undistortedPts(1,:);
  diffVec = undistortedPts( size(undistortedPts,1), : ) - line.pt;
  line.vec = diffVec / norm( diffVec );

  %midPt = line.pt + 0.5 * diffVec;
  %diffFromMid = bsxfun( @minus, undistortedPts, midPt );
  %dists2Pts = norms( diffFromMid, 2, 2 );
  %[~,closestPt2MidIndx] = min( dists2Pts );
  %closestPt2Mid = undistortedPts(closestPt2MidIndx,:);
  %dists = findDistsBetweenPtsAndLine( closestPt2Mid', line );
  %err = norm( dists ) / numel( dists );
  %err = dist ./ norm( diffVec );

  dists = findDistsBetweenPtsAndLine( undistortedPts', line );
  err = max( dists ) / norm( diffVec );
  %err = norm( dists(2:end-1) ) / norm( diffVec ) / ( numel(dists)-2 );
  %err = norm( dists(2:end-1) ) / ( numel(dists)-2 );

end
