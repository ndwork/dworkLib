
function [ outLine, matchIndx ] = findLineInFile( file, exp )
  % [ outLine, matchIndx ] = findLineInFile( file, exp )
  %
  % Finds the first line in the file that matches expression
  %
  % Inputs:
  % file - a string with the file to open
  % exp - the experession to match
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

  fp = fopen( file, 'r' );

  outLine = [];
  while feof( fp ) == 0
    thisLine = fgetl( fp );
    matchIndx = regexp( thisLine, exp );
    if matchIndx
      outLine = thisLine;
      fclose( fp );
      break;
    end
  end

  fclose( fp );
end
