
function [ out, u ] = proxNucNorm( in, thresh, varargin )
  % [ out, u ] = proxNucNorm( in, thresh )
  %
  % Returns the proximal operator of the nuclear norm of the input matrix
  %
  % Inputs:
  % in - an input matrix
  % thresh - the thresholding value
  %
  % Optional Inputs:
  % u - the u matrix of the SVD decomposition to make code faster for tall and skinny
  %     matrices
  %
  % Outputs:
  % out - vector
  % u - the u matrix of the svd
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if thresh == 0
    out = in;
    return;
  end

  p = inputParser;
  p.addParameter( 'u', [], @isnumeric );
  p.parse( varargin{:} );
  u = p.Results.u;

  if numel( u ) > 0
    % Find s and v from svd of in*in and use previous u
    % This should be faster for tall and skinny matrices
    [~,s,v] = svd( in' * in, 'econ', 'vector' );
    s = sqrt( s );
  else
    [u,s,v] = svd( in, 'econ', 'vector' );
  end

  s = softThresh( s, thresh );  % Singular values are always real

  out = u * diag(s) * v';

end
