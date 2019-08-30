
function M = mri_ab2Matrix( a, b )
  % Converts the Cayley-Klein parameters into a matrix so that
  % next([ Mxy; conj(Mxy); Mz]) = Matrix * [ Mxy; conj(Mxy); Mz ]
  %
  % Inputs
  % a - a 1D array of the alpha CK parameters
  % b - a 1D array of the beta CK parameters
  %
  % Outputs:
  % M - a 3D array of size 3x3xN where N is the number of elements in a
  %
  % Written by Nicholas Dwork - Copyright 2018
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if numel(a) ~= numel(b), error('a and b must be the same size'); end

  ca = conj(a);
  cb = conj(b);
  
  M = zeros(3,3,numel(a));

  M(1,1,:) = ca.*ca;    M(1,2,:) = -b.*b;    M(1,3,:) = 2*ca.*b;
  M(2,1,:) = -cb.*cb;   M(2,2,:) = a.^2;     M(2,3,:) = 2*a.*cb;
  M(3,1,:) = -ca.*cb;   M(3,2,:) = -a.*b;    M(3,3,:) = a.*ca-b.*cb;
end
