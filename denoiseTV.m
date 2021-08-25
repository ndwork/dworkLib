
function [out,oValues] = denoiseTV( in, lambda, varargin )
  % out = denoiseTV( in, lambda [, varargin ] )
  %
  % perform total variation denoising on the input
  %
  % Inputs:
  % in - input array to be denoised
  % lambda - the regularization parameter
  %
  % Optional Inputs:
  % varargin - the same optional parameters available to lsqrTV
  %
  % Outputs:
  % out - the denoised array
  %
  % Written by Nicholas Dwork - Copyright 2021
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.
  
  applyA = @(x,op) x;

  if nargout > 1
    [out,oValues] = lsqrTV( applyA, in, in, lambda, varargin{:} );
  else
    out = lsqrTV( applyA, in, in, lambda, varargin{:} );
  end
end
