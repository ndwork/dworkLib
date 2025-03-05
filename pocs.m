
function out = pocs( in, projections, varargin )
  % out = pocs( in, projections [, nMaxIter ] )
  %
  % Projection Onto Convex Sets
  %
  % Inputs:
  % in - the input size
  % projections - a cell array of function handles to projection operations
  %
  % Optional Inputs:
  % 
  %
  % Written by Nicholas Dwork, Copyright 2025
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addOptional( 'nMaxIter', 100, @ispositive );
  p.parse( varargin{:} );
  nMaxIter = p.Results.nMaxIter;

  out = in;
  for iter = 1 : nMaxIter
    for pIndx = 1 : numel( projections )
      out = projections{ pIndx }( out );
    end
  end

end
