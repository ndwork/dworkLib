
function A = makeMatrixForLinearTrans( x, trans, varargin )
  % A = makeMatrixForLinearTrans( x, trans [, parallel ] )
  %
  % Inputs:
  % x - an input for the transformation (has the right size; values don't
  %     matter )
  % trans - a function handle for the transformation
  %
  % Optional Inputs:
  % parallel - if set, does parallel processing
  %
  % Outputs:
  % A - the matrix representing the linear transformation
  %
  % Written by Nicholas Dwork - Copyright 2017
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addOptional( 'parallel', 0, @isnumeric );
  p.parse( varargin{:} );
  parallel = p.Results.parallel;

  N = numel(x);
  M = numel( trans(x) );
  A = zeros(M,N);

  if parallel == 0
    for i=1:N
      tmpIn = zeros(size(x));
      tmpIn(i) = 1;
      tmpOut = trans( tmpIn );
      A(:,i) = tmpOut(:);
    end
  else
    parfor i=1:N
      tmpIn = zeros(size(x));
      tmpIn(i) = 1;
      tmpOut = trans( tmpIn );                                                             %#ok<PFBNS>
      A(:,i) = tmpOut(:);
    end
  end

end
