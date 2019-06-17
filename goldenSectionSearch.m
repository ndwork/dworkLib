
function out = goldenSectionSearch( f, LB, UB, varargin )
  % out = goldenSectionSearch( f, LB, UB [, 'tol', tol, 'nMax', nMax ] )
  %
  % finds the minimal point of the function f using a binary search
  % Written according to the notes written by Wotao Yin at
  % http://www.math.ucla.edu/~wotaoyin/math273a/slides/Lec3a_1d_search_273a_2015_f.pdf
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

  R = 0.61803398874989484;  % golden ratio

  a0 = LB;  b0 = UB;
  D = R * ( b0 - a0 );
  a1 = b0 - D;  fa = f( a1 );
  b1 = a0 + D;  fb = f( b1 );

  iterationIndx = 1;  % 1 iteration already done above
  while iterationIndx < nMax
    if b0 - a0 < tol, break; end

    if fa <= fb
      % the minimal point is in [a0, b1]
      b0 = b1;
      b1 = a1;
      fb = fa;
      a1 = b0 - R * ( b0 - a0 );
      fa = f( a1 );
    else
      % the minimal point is in [a1, b0]
      a0 = a1;
      a1 = b1;
      fa = fb;
      b1 = a0 + R * ( b0 - a0 );
      fb = f( b1 );
    end

    iterationIndx = iterationIndx + 1;
  end

  mid = 0.5 * ( a0 + b0 );
  out = mid;
end

