
function displayNice( dispFunctionH, in1, varargin )
  % disFunctionH - the function handle of the display function
  %   (may be plot or semilogy)
  %
  % This function is called by plotnice and semilogynice
  %
  % Written by Nicholas Dwork, Copyright 2019
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

  if nargin-1 - 2*numel(ax) > 2

    if isnumeric( varargin{1} ) || islogical( varargin{1} )

      if mod( numel(varargin)-1, 2 ) == 0
        if numel( ax ) > 0
          dispFunctionH( ax, in1, varargin{1}, 'LineWidth', 1.5, varargin{2:end} );
        else
          dispFunctionH( in1, varargin{1}, 'LineWidth', 1.5, varargin{2:end} );
        end
      else
        if numel( ax ) > 0
          dispFunctionH( ax, in1, varargin{1}, varargin{2}, 'LineWidth', 1.5, varargin{3:end} );
        else
          dispFunctionH( in1, varargin{1}, varargin{2}, 'LineWidth', 1.5, varargin{3:end} );
        end
      end

    else

      if mod( numel(varargin), 2 ) == 0
        if numel( ax ) > 0
          dispFunctionH( ax, in1, 'LineWidth', 1.5, varargin{:} );
        else
          dispFunctionH( in1, 'LineWidth', 1.5, varargin{:} );
        end
      else
        if numel( ax ) > 0
          dispFunctionH( ax, in1, varargin{1}, 'LineWidth', 1.5, varargin{2:end} );
        else
          dispFunctionH( in1, varargin{1}, 'LineWidth', 1.5, varargin{2:end} );
        end
      end

    end

  elseif nargin-1 - 2*numel(ax) > 1

    if numel( ax ) > 0
      dispFunctionH( ax, in1, varargin{1}, 'LineWidth', 1.5 );
    else
      dispFunctionH( in1, varargin{1}, 'LineWidth', 1.5 );
    end

  else

    if numel( ax ) > 0
      dispFunctionH( ax, in1, 'LineWidth', 1.5 );
    else
      dispFunctionH( in1, 'LineWidth', 1.5 );
    end

  end

  if numel( ax ) > 0
    set( ax, 'fontsize', 14, 'LineWidth', 1.5 );
  else
    set( gca, 'fontsize', 14, 'LineWidth', 1.5 );
  end

  addToolbarExplorationButtons( gcf );  % Restore the missing toolbar
end
