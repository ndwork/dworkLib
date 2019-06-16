

function out = goldenSectionSearch( f, LB, UB, varargin )
  % out = goldenSectionSearch( f, LB, UB [, 'tol', tol, 'nMax', nMax ] )
  %
  % finds the minimal point of the function f using a binary search
  %
  % Inputs:
  % f - function handle
  % LB - the lower bound of the root
  % UB - the upper bound of the root
  %
  % Optional Inputs:
  % tol - the tolerance to use when finding the root (default is 1d-6)
  % nMax - the maximum number of iterations (default is 1000)
  %   Note: nMax can equal Inf
  %
  % Written by Nicholas Dwork, Copyright 2019
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addParameter( 'nMax', 1000, @ispositive );
  p.addParameter( 'tol', 1d-6, @ispositive );
  p.parse( varargin{:} );
  nMax = p.Results.nMax;
  tol = p.Results.tol;

  R = 1.61803398874989484;  % golden ratio

  xL = LB;  fL = f( thisLB );
  xU = UB;  fU = f( thisUB );

  D = R * ( xU - xL );
  x1 = xU - D;  f1 = f( x1 );
  x2 = xL + D;  f2 = f( x2 );

  i = 1;  % 1 iteration already done above
  while i < nMax
    if xU - xL < tol, break; end

    if f1 < f2
      xU = x2;
      x2 = x1;
      f2 = f1;
      x1 = xU - R * ( xU - xL )
      f1 = f( x1 );
    else
      xL = x1;
      x1 = x2;
      f1 = f2;
      x2 = xL + R * ( xU - xL );
      f2 = f( x2 );
    end

    i = i + 1;
  end

  mid = 0.5 * ( xL + xB );
  out = mid;
end

