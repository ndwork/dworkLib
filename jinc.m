
function out = jinc(r)
  % out = jinc( r )
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  out = besselj(1,pi*r) ./ (2*r); 
  out(isnan(out)) = 0.7854;
end


