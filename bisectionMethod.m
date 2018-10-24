
function biRoot = bisectionMethod( biFunc, xLB, xUB, varargin )
  % biRoot = bisectionMethod( biFunc, x1, x2, varargin )
  %
  % Finds the root of the function biFunc
  %
  % Inputs:
  % biFunc - function handle; root will be found of this function
  % xLB - lower bound on root
  % xUB - upper bound on root
  %
  % Optional Inputs:
  % accuracy - error in estimate will be less than this amount (or nMax
  %   will have been reached)
  % nMax - the maximum number of iterations to perform
  %
  % Outputs:
  % biRoot - the estimate of the root
  %
  % Written by Nicholas Dwork - Copyright 2018
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addOptional( 'accuracy', 1d-16, @isnumeric );
  p.addParameter( 'nMax', 1000, @isnumeric );
  p.parse( varargin{:} );
  accuracy = p.Results.accuracy;
  nMax = p.Results.nMax;

  if accuracy < 0, error( 'Accuracy must be greater than 0' ); end;
  if nMax < 1, error('nMax must be greater than 1' ); end;

  f = biFunc( xLB );
  fMid = biFunc( xUB );
  if f*fMid >= 0
    error([ 'The root must be bracketed by the lower', ...
            'bound and the upper bound.' ]);
  end

  if f < 0
    dx = xUB-xLB;  biRoot = xLB;
  else
    dx = xLB-xUB;  biRoot = xUB;
  end

  for n=1:nMax
    dx = dx * 0.5;
    xMid = biRoot + dx;
    fMid = biFunc( xMid );
    if fMid < 0, biRoot = xMid; end;
    if abs(dx) < accuracy || fMid == 0, return, biRoot; end;
	end

end

