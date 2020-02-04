
function out = scaleImg( img, varargin )
  % out = scaleImg( img [, outMinMax, inMinMax ] )
  %
  % Function scaled image affinely so that inMinMax is scaled to outMinMax.
  % If img is complex, it scales the magnitude.
  %
  % Inputs:
  % outMinMax - 2 element array specifying output min and max
  %   or an array from which the min and max are derived.
  %   If not specified, default is [0 1]
  % inMinMax - 2 element array specifying input min and max
  %   or an array from which the min and max are derived
  %   If not specified, true image statistics are used
  %
  % Outputs:
  % out - output image
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It 
  % is offered without any warranty expressed or implied, including the 
  % implied warranties of merchantability or fitness for a particular 
  % purpose.

  if nargin < 1
    disp( 'Usage:  out = scaleImg( img [, outMinMax, inMinMax ] )' );
    return;
  end

  defaultInMinMax = [ min(img(:)) max(img(:)) ];
  defaultOutMinMax = [0 1];

  p = inputParser;
  p.addOptional( 'outMinMax', defaultOutMinMax );
  p.addOptional( 'inMinMax', defaultInMinMax );
  p.parse( varargin{:} );
  inMinMax = p.Results.inMinMax;
  outMinMax = p.Results.outMinMax;

  if numel( inMinMax ) == 0, inMinMax = defaultInMinMax; end
  if numel( outMinMax ) == 0, outMinMax = defaultOutMinMax; end
  
  if isempty( inMinMax ), inMinMax=[ min(img(:)) max(img(:)) ]; end
  inMin = min( inMinMax(:) );  outMin = min( outMinMax(:) );
  inMax = max( inMinMax(:) );  outMax = max( outMinMax(:) );

  if inMax < inMin, error('Max must be larger than min'); end
  if outMax < outMin, error('Max must be larger than min'); end

  if min( isreal( img(:) ) ) == 1  % real data
    out = scaleRealImg( img, outMin, outMax, inMin, inMax );
  else  % complex data
    scaledMag = scaleRealImg( abs( img ), outMin, outMax, inMin, inMax );
    out = scaledMag .* exp( 1i * angle( img ) );
  end

end

function out = scaleRealImg( img, outMin, outMax, inMin, inMax )
  if outMax == outMin
    out = outMax * ones( size( img ) );
  elseif inMin == inMax
    if inMin <= outMin
      out = outMin * ones( size( img ) );
    elseif inMin >= outMax
      out = outMax * ones( size( img ) );
    else
      out = inMin * ones( size( img ) );
    end
  else
    tmp = (img - inMin) / ( inMax - inMin );
    out = tmp * (outMax - outMin) + outMin;
    out = min( max( out, outMin ), outMax );
  end
end

