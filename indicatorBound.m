
function out = indicatorBound( in, bound )
  % out = indicatorBound( in, bound )
  %
  % Implements the indicator function to verify that all elements are within a bound
  % out = indicator( -bound <= in <= bound );
  %
  % Inputs:
  % in - an array of values, real or complex
  % bound - a real scalar value specifying the bound
  %
  % Outupts:
  % out - 0 if all values are within bound and Inf otherwise
  %
  % Written by Nicholas Dwork, Copyright 2021
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if max( abs(in(:)) ) > bound
    out = Inf;
  else
    out = 0;
  end
end
