
function out = scatternice( in1, in2, varargin )
  % scatternice( in1 [, in2, 'ax', ax, options ] )
  %
  % Inputs:
  %   in1 - domain values to plot
  %   in2 - range values to plot
  %
  % Optional Inputs:
  %   ax - the axis to plot onto (used with subplot)
  %   options - all optional arguments that scatter accepts
  %
  % Written by Nicholas - Copyright 2016
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp( 'Usage:  scatternice( in1 [, in2, ''ax'', ax, options ] ) ' );
    if nargout > 0, out = []; end
    return;
  end

  ax = [];
  if numel( varargin ) > 0
    axIndx = find( strcmp( varargin, 'ax' ) );
    if numel( axIndx ) > 0
      ax=varargin{axIndx+1};
      varargin = { varargin{1:axIndx-1} varargin{axIndx+2:end} };
    end
  end

  if nargin < 2
    in2 = in1;
    in1 = 1 : numel( in2 );
  end
  
  if numel( varargin ) > 0
    out = scatter( in1, in2, varargin{:}, 'filled' );
  else
    out = scatter( in1, in2, 'filled' );
  end

  if numel( ax ) > 0
    set( ax, 'fontsize', 14 );
  else
    set( gca, 'fontsize', 14 );
  end

  addToolbarExplorationButtons( gcf );  % Restore the missing toolbar
end
