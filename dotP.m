
function out = dotP( a, b )
  % out = dotP( a, b );
  % Computes the dot product between two vectors a and b.
  % a . b = sum( a .* conj(b) )
  %
  % Written by Nicholas Dwork
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if numel( a ) ~= numel( b )  &&  ( numel( a ) ~= 1 && numel( b ) ~= 1 )
    error( 'Inputs to dot product must have the same number of elements.' );
  end

  out = b(:)' * a(:);
end
