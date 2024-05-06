
function plotnice( in1, varargin )
  % plotnice( in1 [, in2, 'ax', ax, options ] )
  %
  % Inputs:
  %   in1 - 1D array to plot
  %   in2 - if in2 is supplied, in1 are the domain values and in2 are the
  %     range values
  %
  % Optional Inputs:
  %   ax - the axis to plot onto (used with subplot)
  %   options - all optional arguments that plot accepts
  %
  % Written by Nicholas - Copyright 2016
  %
  % www.github.com/ndwork/dworkLib
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp( 'Usage:  plotnice( in1 [, in2, ''ax'', ax, options ] )' );
    return;
  end

  displayNice( @plot, in1, varargin{:} );
  return
end
