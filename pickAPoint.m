
function pt = pickAPoint( varargin )
  % pt = pickAPoint( [ arg1, arg2, 'range', range ] )
  %
  % Find the coordinates of a point in an image
  %
  % Opitonal Inputs:
  % arg1 - can either be an image or a scalar specifying the scale of the
  %   display
  % arg2 - can either be an image or a scalar specifying the scale of the
  %   display
  % 'range' - the dynamic range to use in the display
  %
  % Outputs:
  % pt - [x y]
  %
  % Written by Nicholas Dwork, Copyright 2018
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  isNumericOrLogical = @(x) isnumeric(x) || islogical(x);

  p = inputParser;
  p.addOptional( 'arg1', [], isNumericOrLogical );
  p.addOptional( 'arg2', [], isNumericOrLogical );
  p.addParameter( 'range', [], isNumericOrLogical );
  p.parse( varargin{:} );
  arg1 = p.Results.arg1;
  arg2 = p.Results.arg2;
  range = p.Results.range;

  img = [];
  scale = [];
  if numel( arg1 ) > 1
    img = arg1;
  elseif numel( arg1 ) == 1
    scale = arg1;
  end
  
  if numel( arg2 ) > 1
    if numel( img ) > 0, error('Img already specified'); end
    img = arg2;
  elseif numel( arg2 ) == 1
    if numel( scale ) > 0, error('scale already specified'); end
    scale = arg2;
  end

  if numel( scale ) == 0, scale = 1; end

  if numel( img ) > 0
    figure; imshowscale( img, scale, 'range', range );
  end

  pt = getFeaturesFromImg( 1, scale );
  showFeaturesOnImg( pt, 'scale', scale );
end
