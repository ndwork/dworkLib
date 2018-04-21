
function out = projV1ontoV2( v1, v2 )
  % out = projectV1ontoV2( v1, v2 )
  % Projects vector v1 onto v2
  %
  % Inputs:
  % v1,v2 - vectors (1D arrays) of the same size
  %
  % Outputs:
  % out - a vector the size of v1
  %
  % Written by Nicholas Dwork - Copyright 2018
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 2
    disp('Usage:  out = projectV1ontoV2( v1, v2 )');
    return
  end

  out = dotP( v1, v2 ) / norm( v2, 2 ).^2 .* v2;
end
