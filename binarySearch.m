
function out = binarySearch( f, LB, UB, varargin )
  % out = binarySearch( f, LB, UB [, 'tol', tol, 'nMax', nMax, 'verbose', verbose ] )
  %
  % finds the root of the function f using a binary search
  %
  % Inputs:
  % f - function handle
  % LB - the lower bound of the root
  % UB - the upper bound of the root
  %
  % Optional Inputs:
  % tol - the tolerance to use when finding the root (default is 1d-6)
  % nMax - the maximum number of iterations (default is 1000)
  %
  % Written by Nicholas Dwork, Copyright 2019
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 3
    disp( 'Usage: out = binarySearch( f, LB, UB [, ''tol'', tol, ''nMax'', nMax ] )' );
    return
  end

  p = inputParser;
  p.addRequired( 'LB', @isfinite );
  p.addRequired( 'UB', @isfinite );
  p.addParameter( 'nMax', 1000, @ispositive );
  p.addParameter( 'tol', 1d-6, @ispositive );
  p.addParameter( 'verbose', false, @(x) islogical(x) || isnueric(x) );
  p.parse( LB, UB, varargin{:} );
  nMax = p.Results.nMax;
  tol = p.Results.tol;
  verbose = p.Results.verbose;

  thisLB = LB;  fLB = f( thisLB );
  thisUB = UB;  %fUB = f( thisUB );  % fUB is never used

  i = 0;
  while i < nMax
    if verbose == true
      disp([ 'binarySearch: Working on iteration ', num2str(i), ' of ', num2str(nMax) ' possible.' ]);
    end
    mid = 0.5 * ( thisUB + thisLB );
    if thisUB - thisLB < tol, break; end

    fMid = f( mid );
    if fMid == 0, break; end

    if sign( fLB ) == sign( fMid )
      thisLB = mid;
      fLB = fMid;
    else
      thisUB = mid;
    end

    i = i + 1;
  end

  out = mid;
end

