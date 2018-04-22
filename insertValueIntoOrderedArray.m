
function newArray = insertValueIntoOrderedArray( array, value )
  % newArray = insertValueIntoOrderedArray( array, value )
  %
  % Insert a new value into an ordered array so that the new array remains
  % ordered.
  % If the value already exists in the array, nothing is inserted
  %
  % Inputs:
  % array - a 1D (column or row) array
  % value - the value to insert
  %
  % Written by Nicholas Dwork - Copyright 2017
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if numel( array ) == 0
    newArray = value;
    return
  end

  [~,minIndx] = min( abs( array - value ) );

  if iscolumn( array )

    if array( minIndx ) - value < 0
      newArray = [ array(1:minIndx); 0; array(minIndx+1:end); ];
    elseif array( minIndx ) - value > 0
      newArray = [ array(1:minIndx-1); 0; array(minIndx:end); ];
    else
      newArray = array;
    end

  else  % row array

    if array( minIndx ) - value < 0
      newArray = [ array(1:minIndx) 0 array(minIndx+1:end) ];
    elseif array( minIndx ) - value > 0
      newArray = [ array(1:minIndx-1) 0 array(minIndx:end) ];
    else
      newArray = array;
    end

  end

end
