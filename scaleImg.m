% out = scaleImg( img [, inMinMax, outMinMax ] )
%
% Function scaled image affinely so that inMinMax is scaled to outMinMax
%
% Inputs:
% inMinMax - 2 element array specifying input min and max
%   If not specified, true image statistics are used
% outMinMax - 2 element array specifying output min and max
%   If not specified, default is [0 1]
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

function out = scaleImg( img, varargin )
  defaultInMinMax = [ min(img(:)) max(img(:)) ];
  defaultOutMinMax = [0 1];

  p = inputParser;
  p.addOptional( 'inMinMax', defaultInMinMax );
  p.addOptional( 'outMinMax', defaultOutMinMax );
  p.parse( varargin{:} );
  inMinMax = p.Results.inMinMax;
  outMinMax = p.Results.outMinMax;

  inMin = inMinMax(1);  outMin = outMinMax(1);
  inMax = inMinMax(2);  outMax = outMinMax(2);

  if inMax <= inMin, error('Max must be larger than min'); end;
  if outMax <= outMin, error('Max must be larger than min'); end;

  tmp = (img - inMin) / ( inMax - inMin );
  out = tmp * (outMax - outMin) + outMin;
end
