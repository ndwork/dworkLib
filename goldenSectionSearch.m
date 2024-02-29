
function [out,iterationIndx] = goldenSectionSearch( f, LB, UB, varargin )
  % [out,iterationIndx] = goldenSectionSearch( f, LB, UB [, ...
  %   'tol', tol, 'nMax', nMax, 'verbose', true/false ] )
  %
  % finds the minimal point of the function f using a binary search
  % Written according to the notes written by Wotao Yin at
  % http://www.math.ucla.edu/~wotaoyin/math273a/slides/Lec3a_1d_search_273a_2015_f.pdf
  % Another reference:  https://xdze2.github.io/goldsearch.html
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
  % Outputs:
  % out - result of optimization parameter
  %
  % Optional output:
  % iterationIndx - the iteration index when the optimization ended
  %
  % Written by Nicholas Dwork, Copyright 2019
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  defaultNMax = 1000;
  defaultTol = 1d-6;

  p = inputParser;
  p.addParameter( 'nMax', defaultNMax, @ispositive );
  p.addParameter( 'tol', defaultTol, @ispositive );
  p.addParameter( 'verbose', 0, @(x) ispositive(x) || islogical(x) );
  p.parse( varargin{:} );
  nMax = p.Results.nMax;
  tol = p.Results.tol;
  verbose = p.Results.verbose;

  if numel( nMax ) == 0, nMax = defaultNMax; end
  if numel( tol ) == 0, tol = defaultTol; end

  R = 0.61803398874989484;  % golden ratio

  a0 = LB;  b0 = UB;
  D = R * ( b0 - a0 );
  a1 = b0 - D;  
  b1 = a0 + D; 

  fValues = cell(2,1);
  parfor i = 1 : 2
    if i == 1
      fValues{i} = f( a1 );   %#ok<PFBNS>
    else
      fValues{i} = f( b1 );
    end
  end
  fa = fValues{1};
  fb = fValues{2};
  %fb = f( b1 );
  %fa = f( a1 );

  iterationIndx = 1;  % 1 iteration already done above
  while iterationIndx <= nMax
    if 0.5 * ( b0 - a0 ) < tol, break; end
    if verbose ~= false
      disp([ 'goldenSectionSearch: Working on iteration ', num2str(iterationIndx), ...
        ',  current error: ', num2str( 0.5 * ( b0 - a0 ) ) ]);
    end

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

  iterationIndx = iterationIndx - 1;
  if verbose ~= false
    disp([ 'goldenSectionSearch has completed in ', num2str(iterationIndx), ...
      ' iterations.' ]);
  end

  mid = 0.5 * ( a0 + b0 );
  out = mid;
end

