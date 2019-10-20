

function out = computeGradient( in, varargin )
  % out = computeGradient( in, varargin )
  %
  % This function computes the gradient (or the derivative) of the input
  % with circular boundary conditions
  
  % Inputs:
  % in - a multi-dimensional array
  %
  % Output:
  % out - an array of dimension equal to dimension of input plus one where the gradient
  %   dimension is the last dimension of the output.
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  p = inputParser;
  p.addParameter( 'op', [], @(x) true );
  p.parse( varargin{:} );
  op = p.Results.op;

  if strcmp( 'transp', op )
    sIn = size( in );
    out = zeros( sIn(1:end-1) );

    for dimIndx = 1 : sIn(end)
      cmd2run = [ 'subIn = in(', repmat(':,', [1 ndims(in)-1]), num2str(dimIndx), ');' ];
      eval( cmd2run );
      STin = circshift( subIn, -1, dimIndx );
      out = out + ( STin - subIn );
    end

  else
    out = zeros( [ size( in ) ndims(in) ] );

    for dimIndx = 1 : ndims(in)
      SImg = circshift( in, 1, dimIndx );  % Shifted image
      dimD = SImg - in;   %#ok<NASGU>
      cmd2run = [ 'out(', repmat(':,', [1 ndims(in)]), num2str(dimIndx), ') = dimD;' ];
      eval( cmd2run );
    end

  end
end
