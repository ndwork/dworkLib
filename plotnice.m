
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

  defaultIn2 = [];
  p = inputParser;
  p.addOptional( 'in2', defaultIn2 );
  p.parse( varargin{:} );
  in2 = p.Results.in2;

  if numel(in2)>0
    plot( in1, in2, 'LineWidth', 2, varargin{2:end} );
  else
    plot( in1, 'LineWidth', 2, varargin{:} );
  end
  set( gca, 'fontsize', 18, 'LineWidth', 1.5 );
end
