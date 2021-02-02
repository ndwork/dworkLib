
function out = shiftImg( img, shifts )
  % out = shiftImg( img, shifts )
  %
  % Shift the image; shifts are defined according to circshift
  % Regions with unknown data are zero filled.
  %
  % Inputs:
  % img - an array of any number of dimensions
  % shifts - a 2 element array specifying vertical and horiztonal shift
  %
  % Outputs:
  % out - the shifted image
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  out = circshift( img, shifts );

  for dim = 1 : ndims( img )
    if shifts( dim ) == 0
      continue;

    elseif shifts( dim ) > 0
      cmd = [ 'out( ', repmat( ':, ', [ 1, dim-1 ] ), '1:shifts(', num2str(dim), ')', ...
        repmat( ', :', [ 1, ndims(img)-dim] ), ' ) = 0;' ];
      eval( cmd );

    elseif shifts( dim ) < 0
      cmd = [ 'out( ', repmat( ':, ', [ 1, dim-1 ] ), 'end+shifts(', num2str(dim), ')+1:end', ...
        repmat( ', :', [ 1, ndims(img)-dim] ), ' ) = 0;' ];
      eval( cmd );
    end
  end

end

