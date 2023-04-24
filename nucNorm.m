
function out = nucNorm( in )
  % out = nucNorm( in )
  %
  % Computes the nuclear norm of the input.  The nuclear norm is the input
  % of the singular values of the matrix.
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  [~,s,~] = svd( in, 'econ', 'vector' );
  out = sum( s );
end

