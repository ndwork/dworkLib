
function out = stemnice( in1, varargin )
  % out = stemnice( in1 [, in2, 'ax', ax, options ] )
  %
  % Inputs:
  %   in1 - 1D array to stem plot
  %   in2 - if in2 is supplied, in1 are the domain values and in2 are the
  %     range values
  %
  % Optional Inputs:
  %   ax - the axis to plot onto (used with subplot)
  %   options - all optional arguments that plot accepts
  %
  % Written by Nicholas - Copyright 2018
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp( 'Usage: out = stemnice( in1 [, in2, ''ax'', ax, options ] )' );
    return
  end

  ax = [];
  if numel( varargin ) > 0
    axIndx = find( strcmp( varargin, 'ax' ) );
    if numel( axIndx ) > 0
      ax=varargin{axIndx+1};
      varargin = { varargin{1:axIndx-1} varargin{axIndx+2:end} };
    end
  end

  if nargin - 2*numel(ax) > 2

    if isnumeric( varargin{1} ) || islogical( varargin{1} )

      if mod( numel(varargin)-1, 2 ) == 0
        if numel( ax ) > 0
          out = stem( ax, in1, varargin{1}, 'LineWidth', 1.5, varargin{2:end} );
        else
          out = stem( in1, varargin{1}, 'LineWidth', 1.5, varargin{2:end} );
        end
      else
        if numel( ax ) > 0
          out = stem( ax, in1, varargin{1}, varargin{2}, 'LineWidth', 1.5, varargin{3:end} );
        else
          out = stem( in1, varargin{1}, varargin{2}, 'LineWidth', 1.5, varargin{3:end} );
        end
      end

    else

      if mod( numel(varargin), 2 ) == 0
        if numel( ax ) > 0
          out = stem( ax, in1, 'LineWidth', 1.5, varargin{:} );
        else
          out = stem( in1, 'LineWidth', 1.5, varargin{:} );
        end
      else
        if numel( ax ) > 0
          out = stem( ax, in1, varargin{1}, 'LineWidth', 1.5, varargin{2:end} );
        else
          out = stem( in1, varargin{1}, 'LineWidth', 1.5, varargin{2:end} );
        end
      end

    end

  elseif nargin - 2*numel(ax) > 1

    if numel( ax ) > 0
      out = stem( ax, in1, varargin{1}, 'LineWidth', 1.5 );
    else
      out = stem( in1, varargin{1}, 'LineWidth', 1.5 );
    end

  else

    if numel( ax ) > 0
      out = stem( ax, in1, 'LineWidth', 1.5 );
    else
      out = stem( in1, 'LineWidth', 1.5 );
    end

  end

  if numel( ax ) > 0
    set( ax, 'fontsize', 14, 'LineWidth', 1.5 );
  else
    set( gca, 'fontsize', 14, 'LineWidth', 1.5 );
  end
end
