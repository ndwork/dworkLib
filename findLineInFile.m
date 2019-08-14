
function [ out, matchIndx ] = findLineInFile( file, exp, varargin )
  % [ outLine, matchIndx ] = findLineInFile( file, exp [, 'nLines', nLines ] )
  %
  % Finds the first line in the file that matches expression
  %
  % Inputs:
  % file - a string with the file to open
  % exp - the experession to match
  %
  % Optional Inputs:
  % nLines - the number of lines to acquire
  %
  % Outputs:
  % outLine - the line that matched expression
  %
  % Optional Outputs:
  % matchIndx - the index of the line where exp was matched
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addParameter( 'nLines', 1, @ispositive );
  p.parse( varargin{:} );
  nLines = p.Results.nLines;
  
  fp = fopen( file, 'r' );

  if nLines > 1
    out = cell( nLines, 1 );
  else
    out = [];
  end

  while feof( fp ) == 0
    thisLine = fgetl( fp );
    matchIndx = regexp( thisLine, exp );
    if matchIndx
      if nLines > 1
        out{1} = thisLine;
        for lineIndx = 2:nLines
          if feof( fp ), break; end
          thisLine = fgetl( fp );
          out{ lineIndx } = thisLine;
        end
      else
        out = thisLine;
      end
      break;
    end
  end

  fclose( fp );
end
