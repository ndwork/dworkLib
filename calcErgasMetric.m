
function ergas = calcErgasMetric( color, monochrome, fused, varargin )
  % ergas = calcErgasMetric( color, monochrome, fused, varargin )
  %
  % Inputs:
  % color - an MxNx3 array, where the third dimension represents color
  % monochrome - an MxN array
  % fused - an MxNx3 array representing the fused image
  %
  % Optional Inputs:
  % h - ergas parameter (result is scaled by h/l)
  % l - ergas parameter (result is scaled by h/l)
  %
  % Outputs:
  % ergas - a scalar value
  %
  % Written by Nicholas Dwork - Copyright 2017
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  defaultH = 1;
  defaultL = 1;
  p = inputParser;
  p.addParameter( 'h', defaultH );
  p.addParameter( 'l', defaultL );
  p.parse( varargin{:} );
  h = p.Results.h;
  l = p.Results.l;

  sImg = size( fused );
  nPix = sImg(1)*sImg(2);
  redMSE = norm( color(:,:,1) - fused(:,:,1), 2 )^2 / nPix;
  greenMSE = norm( color(:,:,2) - fused(:,:,2), 2 )^2 / nPix;
  blueMSE = norm( color(:,:,3) - fused(:,:,3), 2 )^2 / nPix;
  monMSE = norm( monochrome - mean(fused,3), 2 )^2 / nPix;

  tmpRed = color(:,:,1);    meanRed = mean( tmpRed(:) );
  tmpGreen = color(:,:,2);  meanGreen = mean( tmpGreen(:) );
  tmpBlue = color(:,:,3);   meanBlue = mean( tmpBlue(:) );
  meanMon = mean( monochrome(:) );

  ergas = 100*h/l * sqrt( 0.25 * ( ...
    ( redMSE / meanRed )^2 + ( greenMSE / meanGreen )^2 + ...
    ( blueMSE / meanBlue )^2 + ( monMSE / meanMon )^2 ) );

end
