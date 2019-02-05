
function nLines = countNumLinesInFile( file )
  % nLines = countNumLinesInFile( file );
  %
  % Inputs:
  % file - string with filename
  %
  % Outputs:
  % nLines - the number of lines in the file
  %
  % Written by Nicholas Dwork - Copyright 2017
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  fid = fopen(file);

  nLines = 0;
  while feof(fid) == 0
    fgetl(fid);
    nLines = nLines + 1;
  end

  fclose(fid);
end
