
function legHandle = legendnice( varargin )
  % legHandle = legendnice( names [, 'ax', ax ] )
  %
  % Makes a legend with large font and puts it in the best place
  %
  % Inputs:
  % names - a cell array that will be the legend
  %
  % Outputs:
  % legHandle - the handle to the legend
  %
  % Optional Inputs:
  % ax - the axis on which to plot the legend
  %
  % Written by Nicholas - Copyright 2017
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if numel( varargin ) == 0, error('Must supply names'); end
  
  ax = [];
  axIndx = find( strcmp( varargin, 'ax' ) );
  if numel( axIndx > 0 )
    ax=varargin{axIndx+1};
    names = { varargin{1:axIndx-1} varargin{axIndx+2:end} };
  else
    names = varargin;
  end

  if numel( ax ) > 0
    legHandle = legend( ax, names{:}, 'Location', 'best' );
  else
    legHandle = legend( names{:}, 'Location', 'best' );
  end

  set( legHandle, 'FontSize', 18 )
end
