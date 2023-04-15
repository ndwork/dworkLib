
function value = findValueBelowFraction( data, fraction, varargin )
  % value = findValueBelowFraction( data, fraction [, 'tol', tol ] )
  %
  % This function finds the value where fraction amount of the data is
  % above this value.
  % For example, if fraction is 0.35, then this function finds the value
  % where approximately 35% of the data is above this value.
  %
  % Inputs:
  % tol - value is within tol of the answer
  %
  % Written by Nicholas Dwork - Copyright 2017
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addParameter( 'tol', 1d-8, @ispositive );
  p.parse( varargin{:} );
  tol = p.Results.tol;

  if nargin < 1
    disp( 'Usage: value = findValueBelowFraction( data, fraction )' );
    return
  end
  
  nData = numel( data );
  lowerValue = min( data(:) );
  upperValue = max( data(:) );

  f = @( value ) fraction - sum( data(:) >= value ) / nData; 
  value = binarySearch( f, lowerValue, upperValue, 'tol', tol, 'nMax', Inf );
end
