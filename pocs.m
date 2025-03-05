
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
  p.addParameter( 'sequential', false );
  p.parse( varargin{:} );
  nMaxIter = p.Results.nMaxIter;
  sequential = p.Results.sequential;

  out = in;

  nProjections = numel( projections );

  if sequential

    for iter = 1 : nMaxIter
      for pIndx = 1 : nProjections
        out = projections{ pIndx }( out );
      end
    end

  else
    sIn = size( in );
    projectionOutputs = cell( 1, nProjections );
    for iter = 1 : nMaxIter
      parfor pIndx = 1 : nProjections
        projOutput = projections{ pIndx }( out );
        projectionOutputs{ pIndx } = projOutput(:);
      end
      out = reshape( mean( cell2mat( projectionOutputs ), 2 ), sIn );
    end
  end

end
