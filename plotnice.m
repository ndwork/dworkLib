
function plotnice( in1, varargin )
  % plotnice( in1 [, in2, options ] )
  %
  % Inputs
  % in1 - 1D array to plot
  % in2 - if in2 is supplied, in1 are the domain values and in2 are the
  %   range values
  % options - all optional arguments that plot accepts
  %
  % Written by Nicholas - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin > 2

    if isnumeric( varargin{1} )

      if mod( numel(varargin)-1, 2 ) == 0
        plot( in1, varargin{1}, 'LineWidth', 2, varargin{2:end} );
      else
        plot( in1, varargin{1}, varargin{2}, 'LineWidth', 2, varargin{3:end} );
      end

    else

      if mod( numel(varargin), 2 ) == 0
        plot( in1, 'LineWidth', 2, varargin{:} );
      else
        plot( in1, varargin{1}, 'LineWidth', 2, varargin{2:end} );
      end

    end

  elseif nargin > 1

    plot( in1, varargin{1}, 'LineWidth', 2 );

  else
    
    plot( in1, 'LineWidth', 2 );

  end

  set( gca, 'fontsize', 16, 'LineWidth', 1.5 );
end
