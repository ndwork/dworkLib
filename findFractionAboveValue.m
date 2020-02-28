
function value = findFractionAboveValue( data, fraction )
  % value = findFractionAboveValue( data, fraction )
  % This function finds the value where fraction amount of the data is
  % above this value.
  % For example, if fraction is 0.35, then this function finds the value
  % where approximately 35% of the data is above this value.
  %
  % Written by Nicholas Dwork - Copyright 2017
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp( 'Usage: value = findFractionAboveValue( data, fraction )' );
    return
  end
  
  nData = numel( data );
  lowerValue = min( data(:) );
  upperValue = max( data(:) );

  fractionAbove = 0;
  lastFractionAbove = fractionAbove - 1.0;

  % Conduct a binary search
  while( lastFractionAbove ~= fractionAbove )
    lastFractionAbove = fractionAbove;

    value = lowerValue + 0.5 * ( upperValue - lowerValue );
    fractionAbove = sum( data(:) > value ) / nData;

    if fractionAbove == fraction
      return;
    elseif fractionAbove > fraction
      lowerValue = value;
    else
      upperValue = value;
    end
  end

end
