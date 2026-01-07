

function out = fixFileseps( in )
  % out = fixFileseps( in )
  % 
  % Replaces all instances of '/' or '\' in a string with the filesep() character
  %
  % Inputs:
  % in - a string specifying a filename
  %
  % Outputs:
  % out - a string specifying the filename using filesep characters
  %
  % Written by Nicholas Dwork, Copyright 2026
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  filenameParts = split( in, '/' );
  if numel( filenameParts ) == 1   %#ok<ISCL>
    out = in;
  else
    out = strjoin( filenameParts, filesep() );
  end

  filenameParts = split( out, '\' );
  if numel( filenameParts ) > 1
    out = strjoin( filenameParts, filesep() );
  end

end
