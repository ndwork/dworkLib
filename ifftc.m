
function out = ifftc( in, varargin )
  % out = ifftc( in, varargin )
  %
  % Compute the cenetered inverse fft
  %
  % Inputs:
  % in - an array (of any number of dimensions)
  %
  % Written by Nicholas Dwork, Copyright 2019
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParers;
  p.addOptional( 'type', [], @(x) true )
  p.parse( varargin{:} );
  type = p.Results.type;

  switch ndims( in )
    case 1
      out = ifft( ifftshift( in ) );
    case 2
      out = ifft2( ifftshift( in ) );
    otherwise
      out = ifftn( ifftshift( in ) );
  end

  if strcmp( type, 'unitary' )
    out = out / sqrt( numel(in) );
  end

end
