
function out = calcMSE( est, true )
  % Calculates the mean squared error of an estimate
  %
  % out = mse( est, true )
  %
  % Inputs:
  % est - the array of estimate values
  % true - an array of the same size as est that represents the correct values
  %
  % Outputs:
  % out - the mean square error value
  %
  % Written by Nicholas Dwork - Copyright 2021
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 2
    disp( 'Usage:  out = mse( est, true )' );
    if nargout > 0, out = []; end
    return;
  end

  out = norm( est(:) - true(:) )^2 / numel( est );
end
