
function [x1,x2] = rootsOfQuadratic( a, b, c )
  % [x1,x2] = rootsOfQuadratic( a, b, c )
  %
  % Finds the roots of the quadratic polynomial a x^2 + b x + c = 0
  % in a numerically stable way.  This method was taken from:
  % https://people.csail.mit.edu/bkph/articles/Quadratics.pdf
  %
  % Written by Nicholas Dwork - Copyright 2017
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  x1 = ( -b - sign(b) * sqrt( b*b - 4*a*c ) ) / ( 2*a );
  x2 = c / ( a * x1 );
end

