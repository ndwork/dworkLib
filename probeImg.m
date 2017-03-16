
function probeImg( varargin )
  % probeImg( [ n ] )
  %
  % Displays the (y,x) coordinates of clicked points on an image, and 
  % (optionally) the corresponding image value.
  %
  % Inputs:
  % n (optional) - the number of points to probe
  %
  % Written by Nicholas Dwork - Copyright 2017
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addOptional( 'n', 1, @(x) isnumeric(x) && mod(x,1)==0 && x>0 );
  p.parse( varargin{:} );
  n = p.Results.n;

  img = getimage(gcf);

  for i=1:n
    [x,y] = ginput(1);
    v = img(y,x);
    disp([ '(y,x) v:   (', num2str(y), ', ', num2str(x) ')  ', num2str(v) ]);
  end

end
