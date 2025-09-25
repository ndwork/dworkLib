
function out = backProject( sino, projAngles, sImg )
  % out = backProject( sino, projAngles, sImg )
  %
  % Inputs:
  % sino - parallel beam sinogram
  %
  % Written by Nicholas Dwork, Copyright 2025
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if size( sino, 1 ) ~= sImg(1)
    error( 'The first dimension of sino must be the number of rows in the image' );
  end

  if size( sino, 2 ) ~= numel( projAngles )
    error( 'The second dimension of sino must be the same as the number of projection angles' );
  end

  out = zeros( sImg );

  for angleIndx = 1 : numel( projAngles )
    projAngle = projAngles( angleIndx ) * 180/pi;

    proj = repmat( sino(:,angleIndx), [ 1 sImg(2) ] );
    proj = imrotate( proj, -projAngle, 'crop' );

    out = out + proj;
  end

end
