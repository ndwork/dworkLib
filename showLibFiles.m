
function showLibFiles( varargin )
  % showLibFiles( [ libName, pattern ] )
  % libName - an optional argument specifying the library name of interest
  % By default, this function shows all libraries.
  % If a library name is specified, then this function shows all the files
  % in that library.
  %
  % Optional Inputs:
  % libName - an optional imput of the library to list
  % pattern - an optional input specying a search pattern that must be
  %   matched for a file to be listed
  %
  % Written by Nicholas - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  defaultPattern = '.*';
  p = inputParser;
  p.addOptional( 'libName', [], @(x) true );
    % must specify validation function above to keep addOptional from
    % thinking that argument is a parameter's name.
  p.addOptional( 'pattern', defaultPattern, @(x) true );
  p.addOptional( 'inFileName', false, @islogical );
  p.parse( varargin{:} );
  libName = p.Results.libName;
  pattern = p.Results.pattern;
  inFileName = p.Results.inFileName;

  if isempty( libName )
    disp('Usage: showLibFiles( [ libName ] )' );
    disp('Libraries: ');
    showLibs;

  else
    myPath = path;
    pathDirs = strsplit( myPath, ':' );

    % Search in filenames
    found = 0;
    for i=1:numel(pathDirs)
      if ~isempty( regexp( pathDirs{i}, libName, 'ONCE' ) )
        files = dir( pathDirs{i} );
        for j=1:numel(files)
          if isempty( regexp(files(j).name, '^\.', 'ONCE' ) ) && ...
            ~isempty( regexpi(files(j).name, pattern, 'ONCE' ) )

            found = 1;
            disp( files(j).name );
          end
        end
      end
    end
    if found == 0
      disp('The pattern does not match with any filenames.');
    end
    if inFileName == true
      return
    end

    foundWithin = 0;
    if ~strcmp( pattern, defaultPattern )
      % Search within files
      for i=1:numel(pathDirs)
        if ~isempty( regexp( pathDirs{i}, libName, 'ONCE' ) )
          files = dir( pathDirs{i} );
          for j=1:numel(files)
            if ~isempty( regexpi(files(j).name, pattern, 'ONCE' ) ), continue; end
            if ~isempty( regexp(files(j).name, '^\.', 'ONCE' ) ), continue; end

            fid = fopen( [ pathDirs{i}, '/', files(j).name ] );
            while feof(fid) == 0
              thisLine = fgetl(fid);
              if regexpi( thisLine, pattern, 'ONCE' )
                if foundWithin == 0
                  if found == 0
                    fprintf('\nBut, it was found within these files: \n');
                  else
                    fprintf('\nAlso, it was found within these files: \n');
                  end
                  foundWithin = 1;
                end
                disp([ '  ', files(j).name ]);
                break
              end
            end
            fclose(fid);
          end
        end
      end
      if found == 0 && foundWithin == 0
        disp('The pattern is not within any filenames.');
      elseif foundWithin == 0
        disp('The pattern is not within any other filenames.');
      end
    end

  end

end

