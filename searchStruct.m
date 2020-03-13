
function matches = searchStruct( in, el2Match, varargin )
  % matches = searchStruct( in, el2Match [, type ] )
  %
  % Inputs:
  % in - the input struct
  % el2Match - the element to match
  %   If el2Match is a string, then it searches the field names.
  %   If el2Match is numeric, then it searches the values
  %
  % Optional Inputs:
  % type - Specifies the aspects of the structure to search.
  %   Can either be 'names' or 'values'.
  %
  % Written by Nicholas - Copyright 2019
  % This file is part of dworkLib located at www.github.com/ndwork/dworkLib
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 2
    disp( 'Usage:  matches = searchStruct( in, el2Match [, type ] )' );
    return
  end

  p = inputParser;
  p.addOptional( 'type', [], @(x) true );
  p.parse( varargin{:} );
  type = p.Results.type;

  if isstruct( el2Match )
    error( 'Cannot search a struct for another struct' );
  end
  
  if numel( type ) == 0
    type = 'names';
    if isnumeric( el2Match ), type = 'values'; end
  end
  
  if strcmp( type, 'names' ) && isnumeric( el2Match )
    error( 'Can only search names with a string regular expression' );
  end

  names = fieldnames( in );

  matches = {};
  matchIndx = 1;

  for nameIndx = 1 : numel ( names )
    if strcmp( type, 'names' )

      if numel( regexp( names{nameIndx}, el2Match ) ) ~= 0
        matches{matchIndx} = names{nameIndx};   %#ok<AGROW>
        matchIndx = matchIndx + 1;
      end

      thisValue = in.( names{nameIndx} );
      if isstruct( thisValue )

        newMatches = searchStruct( thisValue, el2Match, 'names' );
        for newMatchIndx = 1 : numel( newMatches )
          matches{ matchIndx } = [ names{nameIndx}, '.', newMatches{newMatchIndx} ];   %#ok<AGROW>
          matchIndx = matchIndx + 1;
        end
      end

    elseif strcmp( type, 'values' )

      thisValue = in.( names{nameIndx} );

      if isstruct( thisValue )

        newMatches = searchStruct( thisValue, el2Match, 'values' );
        for newMatchIndx = 1 : numel( newMatches )
          matches{ matchIndx } = [ names{nameIndx}, '.', newMatches{newMatchIndx} ];   %#ok<AGROW>
          matchIndx = matchIndx + 1;
        end

      elseif iscell( thisValue )

        for vIndx = 1 : numel( thisValue )
          found = find( thisValue{vIndx} == el2Match, 1 );
          if numel( found ) ~= 0
            matches{matchIndx} = names{nameIndx};   %#ok<AGROW>
            matchIndx = matchIndx+1;
            break;
          end
        end

      else

        theseValues = in.(names{nameIndx});
        if ischar( theseValues ) && ischar( el2Match )
          found = regexp( theseValues, el2Match );
        elseif ischar( theseValues ) && ~ischar( el2Match )
          found = [];
        elseif ~ischar( theseValues ) && ischar( el2Match )
          found = [];
        else
          found = find( theseValues == el2Match, 1 );
        end
        if numel(found) ~= 0
          matches{matchIndx} = names{nameIndx};   %#ok<AGROW>
          matchIndx = matchIndx+1;
        end
      end

    end

  end

  if matchIndx == 1
    matches = [];
  else
    matches = matches( 1 : matchIndx-1 );
  end

end
