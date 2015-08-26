
function showLibFiles( varargin )
  % showLibFiles( [ libName ] )
  % libName - an optional argument specifying the library name of interest
  % By default, this function shows all libraries.
  % If a library name is specified, then this function shows all the files
  % in that library.

  defaultLib = [];
  p = inputParser;
  p.addOptional( 'libName', defaultLib, @(x) true );
    % must specify validation function above to keep addOptional from
    % thinking that argument is a parameter's name.
  p.parse( varargin{:} );
  libName = p.Results.libName;

  if isempty( libName )
    disp('Usage: showLibFiles( [ libName ] )' );
    disp('Libraries: ');
    showLibs;

  else
    myPath = path;
    pathDirs = strsplit( myPath, ':' );

    for i=1:numel(pathDirs)
      if ~isempty( regexp( pathDirs{i}, libName ) )
        files = dir( pathDirs{i} );
        for j=1:numel(files)
          if isempty( regexp(files(j).name, '^\.' ) )
            disp( files(j).name );
          end
        end
      end
    end

  end

end

