

function namePlottedPoints( x, y, names, varargin )
  % namePlottedPoints( x, y, names [, dx, dy, options ] )
  % Names individual points in a plot; this function can be used with plot
  % or scatter.
  %
  % Inputs:
  % x - the domain values that were plotted
  % y - the range values that were plotted
  % names - a cell array of names
  %
  % Optional Inputs:
  %   dx - the horizontal offset of the names (in pixels)
  %   dy - the vertical offset of the names (in pixels)
  %   options - all other options accepted by text
  %
  % Written by Nicholas - Copyright 2018
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addOptional( 'dx', 0, @isnumeric );
  p.addOptional( 'dy', 15, @isnumeric );
  p.parse( varargin{:} );
  dx = p.Results.dx;
  dy = p.Results.dy;

  ax = gca;
  currentAxUnits = ax.Units;
  set( ax, 'Units', 'pixels' );

  axPos = get( ax, 'position' );
  plotWidthPxls = axPos(3);
  plotHeightPxls = axPos(4);

  axXLim = ax.XLim;  axYLim = ax.YLim;
  plotWidthData = axXLim(2) - axXLim(1);
  plotHeightData = axYLim(2) - axYLim(1);

  dxData = dx / plotWidthPxls * plotWidthData;
  dyData = dy / plotHeightPxls * plotHeightData;
  
  t = text( x+dxData, y+dyData, names, 'FontSize', 16 );
  for i=1:numel(t)
    t(i).HorizontalAlignment = 'center';
  end

  set( ax, 'Units', currentAxUnits )
end
