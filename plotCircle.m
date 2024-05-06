
function h = plotCircle( c, r, varargin )
  % h = plotCircle( c, r, varargin )
  %
  % Inputs:
  % c - a two element array specifying the location of the center
  % r - the radius of the circle
  %
  % Optional Inputs:
  % dAngle - the angular increment used when plotting
  % LineWidth - the line width used when plotting
  %
  % Written by Nicholas - Copyright 2016
  %
  % www.github.com/ndwork/dworkLib
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  p = inputParser;
  p.addParameter( 'dAngle', 1/(2*pi), @ispositive );
  p.addParameter( 'LineWidth', 1.5 );
  p.parse( varargin{:} );
  dAngle = p.Results.dAngle;
  LineWidth = p.Results.LineWidth;

  isHoldIsOn = ishold;
  hold on
  angles = 0 : dAngle : 2*pi + dAngle;
  xs = r * cos( angles ) + c(1);
  ys = r * sin( angles ) + c(2);
  h = plot( xs, ys, 'LineWidth', LineWidth );
  if isHoldIsOn == false,  hold off; end
end

