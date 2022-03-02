
function out = sinc( in )
  % out = sinc( in )
  %
  % Computes the cardinal sine of all elements of the input
  %
  % Inputs:
  % in - an array (of any number of dimensions) of scalar values
  %
  % Outputs:
  % out - an array of size equal to in with the results of since applied to
  %       each element
  %
  % Written by Nicholas Dwork - Copyright 2022
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  out = sin( pi * in ) ./ ( pi * in );
  out( in == 0 ) = 1;

end
