
function [ ks, c ] = findRadialDistortionFromLines( lines, varargin )
  % k = findRadialDistortionFromLines( img, lines [, order ] )
  %
  % The radial distortion is according to section 7.4 of Multiple View Geometry, 2nd edition
  % by Hartley and Zisserman.
  % x' = x_c + L(r)( x - x_c )   and   y' = y_c + L(r)( y - y_c )
  %
  % Note: this function assumes the aspect ratio of point coordinates is 1.
  %
  % Inputs:
  % lines - a 1D cell array where each cell contains a set of points belonging to a line
  %   The element of each cell is a NL x 2 array where NL is the number of points for that line
  %
  % Optional Inputs:
  % order - the order of the parameters to identify (default is 1)
  % c0 - the best estimate of the coordinate of the center of the image
  % fixedCenter (true/false) - if true, then c0 is accepted
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
  p.addParameter( 'c0', [], @(x) isnumeric(x) || numel(x) == 0 );
  p.addParameter( 'fixedCenter', false );
  p.parse( varargin{:} );
  order = p.Results.order;
  c0 = p.Results.c0;
  fixedCenter = p.Results.fixedCenter;

  ks0 = zeros( order, 1 );
  if numel( c0 ) == 0, c0 = [0; 0;]; end

  %err = computeLineError( lines{1}, ks0, c0 );

  function err = computeThisError( x )
    if fixedCenter == true
      cThis = c0;
      ksThis = x;
    else
      cThis = x(1:2);
      ksThis = x(3:end);
    end

    err = computeError( lines, ksThis, cThis );
  end

  if fixedCenter == true
    x0 = ks0;
  else
    x0 = [ c0; ks0; ];
  end
  
  options = optimset( 'TolFun', 1e-12, 'MaxFunEvals', 10000 );
  [x,finalErr] = fminsearch( @computeThisError, x0, options );   %#ok<ASGLU>

  if fixedCenter == true
    c = c0;
    ks = x;
  else
    c = x(1:2);
    ks = x(3:end);
  end

end


function err = computeError( lines, ks, c )
  err = 0;
  for i = 1 : numel( lines )
    err = err + computeLineError( lines{i}, ks, c );
  end
end


function err = computeLineError( pts, ks, c )

  undistortedPts = applyRadialDistortion2Pts( pts, ks, c, 'dir', -1 );

  [ pt, vec ] = findBestLineThroughPoints( undistortedPts );
  line.pt = pt;
  line.vec = vec;

  dists = findDistsBetweenPtsAndLine( undistortedPts', line );
  err = sum( dists.^2 );
end
