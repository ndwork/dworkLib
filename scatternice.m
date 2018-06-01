
function scatternice( in1, in2, varargin )
  % scatternice( in1, in2  [, 'ax', ax, options ] )
  %
  % Inputs:
  %   in1 - domain values to plot
  %   in2 - range values to plot
  %
  % Optional Inputs:
  %   ax - the axis to plot onto (used with subplot)
  %   options - all optional arguments that plot accepts
  %
  % Written by Nicholas - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.


  ax = [];
  if numel( varargin ) > 0
    axIndx = find( strcmp( varargin, 'ax' ) );
    if numel( axIndx > 0 )
      ax=varargin{axIndx+1};
      varargin = { varargin{1:axIndx-1} varargin{axIndx+2:end} };
    end
  end

  if numel( varargin ) > 0
    scatter( in1, in2, 'filled', 'LineWidth', 1.5, varargin{:} );
  else
    scatter( in1, in2, 'filled', 'LineWidth', 1.5 );
  end

  if numel( ax ) > 0
    set( ax, 'fontsize', 14, 'LineWidth', 1.5 );
  else
    set( gca, 'fontsize', 14, 'LineWidth', 1.5 );
  end
end