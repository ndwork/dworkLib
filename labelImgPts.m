
function labelImgPts( pts, varargin )
  % labelImgPts( pts [, 'scale', scale, 'inlierIndxs', inlierIndxs] )
  % This function puts the index of the point onto the image
  % It assumes that the image is already displayed and the figure
  %   is set to the current figure.
  % pts - an Nx2 array where the first column is the x location
  %   and the second column is the y location
  % scale (optional) - scales the points to display correctly on a scaled
  %   image.  Its default value is 1.
  %
  % Written by Nicholas Dwork (c) 2015
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  defaultScale = 1.0;
  defaultInlierIndxs = [];
  p = inputParser;
  p.addParameter('scale',defaultScale);
  p.addParameter('inlierIndxs',defaultInlierIndxs);
  p.parse( varargin{:} );
  scale = p.Results.scale;
  inlierIndxs = p.Results.inlierIndxs;

  scaledPts = round( pts * scale );

  [nPts,~] = size(scaledPts);
  for i=1:nPts
    x = scaledPts(i,1);
    y = scaledPts(i,2);

    text( x-2, y-2, num2str(i), 'color', 'blue', 'FontSize', 20 );

    if numel(inlierIndxs) > 0
      if find( inlierIndxs == i, 1, 'first' )
        text( x, y, num2str(i), 'color', 'yellow', 'FontSize', 20 );
      else
        text( x, y, num2str(i), 'color', 'red', 'FontSize', 20 );
      end
    else
      text( x, y, num2str(i), 'color', 'yellow', 'FontSize', 20 );
    end
  end

  drawnow;
end
