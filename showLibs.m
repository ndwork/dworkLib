
function showLibs( )
  % Shows all libraries available

  myPath = path;
  pathDirs = strsplit( myPath, ':' );

  for i=1:numel(pathDirs)
    if pathDirs{i} == '.'
      continue;
    elseif ~isempty( regexp( pathDirs{i}, '^/Applications/MATLAB_', 'ONCE' ) )
      continue
    else
      disp(pathDirs{i});
    end
  end

end
