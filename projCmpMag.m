
function out = projCmpMag( x, t )
  % out = projCmpMag( x, t )
  %
  % Project x so that its magnitude is less than or equal to t
  %
  % Inputs:
  % x - an array of values to project onto the set {|x_i| <= t for all i}
  % t - scalar >= 0
  %
  % out - array of size(x) with projected magnitudes
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if t < 0
    error( 't must be positive' );

  elseif t == 0
    out = x;  out(:) = 0;

  else

    xMag = abs( x );
    scalingFactors = ones( size(x) );
    scalingFactors( xMag > t ) = 1 ./ xMag( xMag > t );

    out = x .* scalingFactors;

  end
end
