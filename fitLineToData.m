
function [m,b] = fitLineToData( x, y )
  % [m,b] = fitLineToData( x, y ) or
  % [m,b] = fitLineToData( y )
  %
  % This function finds m and b so that ||y-(mx+b)||_2 is minimized
  %
  % Inputs:
  % x - (optional) domain values
  % y - range values
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 2
    y = x;
    x = 1:numel(y);
  end

  A = ones( numel(x), 2 );
  A(:,1) = x(:);
  v = A \ y(:);
  m=v(1);
  b=v(2);
end